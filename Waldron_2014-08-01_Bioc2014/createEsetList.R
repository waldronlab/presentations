library(genefilter)
library(survival)
library(logging)

## -----------------------------------------------------------------------------
## Logging setup
## -----------------------------------------------------------------------------
basicConfig()

## -----------------------------------------------------------------------------
loginfo("Inside script createEsetList.R - inputArgs =")

library(package.name, character.only=TRUE)

loginfo(paste("Loading", package.name, sessionInfo()$otherPkgs[[package.name]]$Version))

## -----------------------------------------------------------------------------
## needed functions
## -----------------------------------------------------------------------------
filterQuantile <- function(object, q){
    if (!identical(q >=0 && q < 1, TRUE))
        stop("require 0 <= q < 1")
    if (!identical(class(object) == "ExpressionSet", TRUE))
        stop("object must be an ExpressionSet")
    gene.sd <- esApply(object,1,sd, na.rm=TRUE)
    gene.quantile <- quantile(gene.sd, probs=q)
    actual.makescutoff <- sum(gene.sd < gene.quantile) / length(gene.sd)
    ##make sure the correct number of genes are getting filtered:
    if (abs(q - actual.makescutoff) > 0.01){
        stop("Not scaling this object, likely pre-scaled.")
    }else{
        object <- object[gene.sd > gene.quantile, ]
    }
    return(object)
}
##recursive intersect function
intersectMany <- function(lst){
    ## Find the intersection of multiple vectors stored as elements of a
    ## list, through a tail-recursive function.
    if (length(lst)==2){
        return(intersect(lst[[1]],lst[[2]]))
    }else{
        return(intersectMany(c(list(intersect(lst[[1]],lst[[2]])),lst[-1:-2])))
    }
}

## -----------------------------------------------------------------------------
##load the esets
## -----------------------------------------------------------------------------
data(list=data(package=package.name)[[3]][,3])

strEsets <- ls(pattern="^.*_eset$")
esets <- list()

## -----------------------------------------------------------------------------
##Explicit removal of datasets:
## -----------------------------------------------------------------------------
if(exists("remove.datasets") && any(strEsets %in% remove.datasets)){
    remove.datasets <- remove.datasets[remove.datasets %in% strEsets]
    loginfo(paste("Removing ", paste(remove.datasets, collapse=", "), " (remove.datasets)"))
    strEsets <- strEsets[!strEsets %in% remove.datasets]
}

## -----------------------------------------------------------------------------
##Explicit removal of samples from specified datasets:
## -----------------------------------------------------------------------------
delim <- ":"   ##This is the delimiter used to specify dataset:sample,
               ##e.g. TCGA_eset:TCGA.24.1927
if(exists("remove.samples")){
    ##split into those with the delimiter and those without:
    remove.samples.delim <- grep(delim, remove.samples, fixed=TRUE, value=TRUE)
    ##over-write remove.samples for later, keeping this variable in
    ##its historical use:
    remove.samples <- grep(delim, remove.samples, fixed=TRUE, invert=TRUE, value=TRUE)
    if(length(remove.samples.delim) > 0){
        datasets <- gsub(paste(delim, ".+", sep=""), "", remove.samples.delim)
        samples <- gsub(paste(".+", delim, sep=""), "", remove.samples.delim)
        remove.samples.delim <- lapply(unique(datasets), function(ds){
            samples[datasets %in% ds]
        })
        names(remove.samples.delim) <- unique(datasets)
    }
}

loginfo("Clean up the esets.")
for (strEset in strEsets){
    eset <- get(strEset)
    ##Remove genes which had a single probe mapping to multiple genes:
    eset <- eset[!grepl("///",featureNames(eset),fixed=TRUE),]

    # Run ComBat
    if (exists("combat") && combat) {
        # workaround bug #12
        if (identical(pubMedIds(eset), "21720365")) {
            eset$batch[match("TCGA.23.1023", sampleNames(eset))] <- 12
        }
        if (sum(is.na(eset$batch)) == 0) {
            tmp <- try(sva::ComBat(exprs(eset),
                mod=model.matrix(~rep(1, ncol(eset))), batch=eset$batch), silent=TRUE)
            if (class(tmp)=="matrix") {
                loginfo(paste("Making ComBat correction to", strEset))
                exprs(eset) <- tmp
            }
        } 
    }

    ##samples to be removed
    remove <- rep(FALSE, ncol(eset))
    ##remove samples without required metadata
    if(exists("meta.required") && length(meta.required) > 0){
        for (varname in meta.required){
            if (varname %in% colnames(pData(eset))){
                remove[ is.na(eset[[varname]]) ] <- TRUE
            }
        }
    }
    ##remove samples not matching regexes
    all.rules <- ls(pattern="rule\\.[0-9]+")
    for (one.rule in all.rules){
        this.remove <- !grepl(get(one.rule)[2], eset[[ get(one.rule)[1] ]])
        if(!strict.checking)
            this.remove[ is.na(eset[[ get(one.rule)[1] ]]) ] <- FALSE
        remove[this.remove] <- TRUE
    }
    ##remove samples pre-specified for removal, that have a dataset specified:
    if(exists("remove.samples.delim")){
        if (strEset %in% names(remove.samples.delim)){
            remove[sampleNames(eset) %in% remove.samples.delim[[strEset]]] <- TRUE
        }
    }
    ##remove samples pre-specified for removal, that did *not* have a dataset specified:
    if(exists("remove.samples"))
        remove[sampleNames(eset) %in% remove.samples] <- TRUE
    ##do the actual removal
    eset <- eset[, !remove]
    if (exists("considered.datasets") && !(strEset %in% considered.datasets))
    {
        loginfo(paste("excluding",strEset,
            "(considered.datasets)"))
        next
    }
    ##include study if it has enough samples and events:
    if (exists("min.number.of.events") && !is.na(min.number.of.events)
        && exists("min.sample.size") && !is.na(min.sample.size)
        && min.number.of.events > 0
        && sum(eset$vital_status == "deceased") < min.number.of.events
        || ncol(eset) < min.sample.size)
    {
        loginfo(paste("excluding",strEset,
            "(min.number.of.events or min.sample.size)"))
        next
    }
    if (exists("min.number.of.genes") && nrow(eset) < min.number.of.genes) {
        loginfo(paste("excluding",strEset,"(min.number.of.genes)"))
        next
    }
    ##filter genes with standard deviation below the required quantile
    if(exists("quantile.cutoff") && quantile.cutoff > 0 && quantile.cutoff < 1){
        eset <- filterQuantile(eset, q=quantile.cutoff)
    }
    ##rescale to z-scores
    if(exists("rescale") && rescale){
        exprs(eset) <- t(scale(t(exprs(eset))))
    }
    loginfo(paste("including",strEset))
##    featureNames(eset) <- make.names(featureNames(eset))  ##should not do this, it is irreversible.
    esets[[strEset]] <- eset
    rm(eset)
}

##optionally take the intersection of genes common to all platforms:
if(exists("keep.common.only") && keep.common.only){
    features.per.dataset <- lapply(esets, featureNames)
    intersect.genes <- intersectMany(features.per.dataset)
    esets <- lapply(esets, function(eset){
        eset <- eset[intersect.genes, ]
        return(eset)
    })
}

ids.with.missing.data <- which(sapply(esets, function(X)
                                      sum(!complete.cases(exprs(X))) > 0))
loginfo(paste("Ids with missing data:", paste(names(ids.with.missing.data),
                                              collapse=", ")))

if (length(ids.with.missing.data) > 0 && exists("impute.missing") && impute.missing) {
    for (i in ids.with.missing.data) {
        require(impute)
        exprs(esets[[i]]) = impute.knn(exprs(esets[[i]]))$data
    }
}


if (exists("add.surv.y") && is.function(add.surv.y)) {
    for (i in 1:length(esets)) {
        esets[[i]]$y = add.surv.y(esets[[i]])
    }
} 
