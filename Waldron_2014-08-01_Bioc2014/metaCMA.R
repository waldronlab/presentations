library(metafor)
library(survHD)
library(survcomp)
library(rmeta)

.combineEsets <- function(esets, y="y",
    probesets=featureNames(esets[[1]]),ComBat=FALSE) {
    if (length(esets)==0) return(NULL)
    X = do.call("cbind", lapply(esets, exprs))
    if (!is.null(y) && class(esets[[1]][[y]])=="Surv") {
        y <- do.call("rbind", lapply(esets, function(x) x[[y]]))
        y <- Surv(y[,1], y[,2])
    } else {
        y <- do.call("c", lapply(esets, function(x) as.character(x[[y]])))
    }
    batch =  as.factor(do.call(c,lapply(1:length(esets), function(i)
        rep(names(esets)[i], ncol(esets[[i]])))))

    if (ComBat) {
        X <- sva::ComBat(X[probesets,], mod = model.matrix(~(rep(1, 
                        length(batch)))), batch = batch)
    }

    l = list(X=t(X[probesets,]),y=y,batch=batch)
    colnames(l$X) = make.names(colnames(l$X), unique = FALSE)
    l
}

.calcRMA <- function(idx, coefs, rma.method="FE") {
    res <- lapply(1:nrow(coefs$c),function(i) try(metafor::rma(yi = coefs$c[i,idx], sei =
        coefs$se[i,idx], method=rma.method)))
    # in case something went wrong, exclude the gene    
     res[which(sapply(res, function(x) class(x)[1]) != "rma.uni")] <-
     list(list(pval=1, b=0))

    res    
}

.createXy <-function(idx, esets, y="y", coefs, n,
rma.method="FE", filter.fun=.defaultFilter, modeltype="compoundcovariate"){ 
    probesets = featureNames(esets[[1]])
    coefficients = NULL
    model = NULL
    pvalues <- NULL

    idx.not.want = -which(sapply(esets,filter.fun))
    # idx is minus validation set
    idx.t = c(idx, idx.not.want)
    # no validation data at all? then use all 
    if (length(idx.t) == 0) idx.t <- 1:length(esets)
    if (length(esets[idx.t]) == 0) {
        warning("No training data passed the filter.")
        return(NULL)
    }
    idx.v = c()
    if (length(idx) > 0) idx.v = -idx
    res.rma = c() 
    if (!is.null(n)) {
        res.rma = .calcRMA(idx.t, coefs, rma.method)
        pvalues = as.numeric(sapply(res.rma, function(r) try(r$pval)))
        #ids = apply(coefs$c[,idx],1,function(x) sum(sign(x))==length(x) || sum(sign(x))==-length(x) )
        names(pvalues) = rownames(coefs$c)
        coefficients = sapply(res.rma, function(x) x$b)
        #pvalues[!ids] = 1
        names(coefficients) = make.names(rownames(coefs$c))
        probesets <- head(order(pvalues),n)
        res.rma <- head(res.rma[order(pvalues)],n)
        names(res.rma) <- names(coefficients[probesets])
        model = new("linearriskscore",
        coefficients=coefficients[probesets],modeltype=modeltype)
    } else {
        cat("No RMA because no n")
    }
    xy = list(train    = .combineEsets(esets[idx.t], y, probesets),
              validate = .combineEsets(esets[idx.v], y, probesets),
              model    = model, idx.t = idx.t, idx.v = idx.v, 
              pvalues  = pvalues,
              rma      = res.rma
    )
    xy
}

metaCMA.train <- function(i, esets, y="y", coefs, n, method=NULL, rma.method="FE",
verbose=FALSE, filter.fun=.defaultFilter, modeltype="compoundcovariate", ...) {
    if (verbose) cat("Dataset", i, "of", length(esets),"...\n")
    tmp <- .createXy(-i,  esets, y, coefs, n, rma.method,
    filter.fun,modeltype=modeltype)
    if (is.null(tmp)) {
        return(list)
    } else if (is.null(method)) {
        fit = tmp$model
    } else {
        fit = do.call(method,args=list(X=tmp$train$X, y=tmp$train$y, ... ))
    }
    list(fit = fit, risk = predict(fit, tmp$validate$X, type="lp"))
}

.defaultCoef <- function(x, y, X,...) {
    if (is.Surv(y)) { 
        fit = coxph(y~x) 
        return(summary(fit)$coefficients[c(1,3)])
    }    
    fit = glm(y~x,...)
    cf = summary(fit)$coefficients
    if (nrow(cf)==2) return(cf[2,1:2])
    warning(paste("Coefficients not found",cf))
    c(NA,NA)
}

metaCMA.coefs <- function(esets, y="y", coef.fun=.defaultCoef, ...) {
    res = lapply(esets, function(X) apply(exprs(X), 1, coef.fun, X[[y]], X, ...))
    coefs = list(c = mapply(cbind, lapply(res, function(x) x[1,])),
                 se = mapply(cbind, lapply(res, function(x) x[2,])))
    rownames(coefs$c)  = featureNames(esets[[1]])
    rownames(coefs$se) = featureNames(esets[[1]])
    # only probe sets without missing values and non-zero se in all datasets
    idx = complete.cases(coefs$c) & apply(coefs$se,1,function(x) sum(x==0)==0)
    coefs$c  = coefs$c[idx,]
    coefs$se = coefs$se[idx,]
    coefs
}

.defaultFilter <- function(eset) {
    ncol(eset) < 75
}

metaCMA.opt <- function(esets, ...) {
    .createXy(esets,idx=c(),...)
}

metaCMA <- function(esets, y="y",  coefs=NULL, n, method=NULL,
rma.method="FE", filter.fun=.defaultFilter, modeltype="compoundcovariate",
... ) {
    if (is.null(coefs)) coefs = metaCMA.coefs(esets, y)
    fits  = lapply(1:length(esets), metaCMA.train, esets=esets, y=y, coefs=coefs,
    n=n, method=method,
    rma.method=rma.method,filter.fun=filter.fun, modeltype=modeltype, ...)
    names(fits) = names(esets)
    list(fits=fits, y=y, rma.method=rma.method)
}

metaCMA.eval <- function(ids, esets, y="y", object, tau=365.25*4,...) {
    lapply(ids, function(i)
    evaluate(object$fits[[i]]$risk, measure=new("UnoC"),
        newy=esets[[i]][[y]], add=list(tau=tau) ))
}

.camelCase <- function(x) paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)), sep="")

.forestplot <- function(esets, y, label, rma.method="FE", measure="hr",
at = NULL, concordance=TRUE, xlab=ifelse(concordance, "Concordance", "Hazard ratio"),...) {
    if (is.null(names(esets))) names(esets) = paste("Study", 1:length(esets))
    names(esets) = gsub("_", " ", names(esets))
    if (concordance) {
        coefs = sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~esets[[i]][[label]]))$concordance)   
        res.rma = metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        if (is.null(at)) at = seq(0.4,1.0,0.1)
        forest.rma(res.rma, 
        xlab=xlab,  slab=sapply(names(esets), .camelCase),
        refline=0.5, at=at,...)
        return(res.rma)
    } else {
        coefs = sapply(1:length(esets), function(i)
        .defaultCoef(esets[[i]][[label]], esets[[i]][[y]],family="binomial"))
        res.rma = metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        if (is.null(at)) at = log(c(0.25,1,4,20))
        forest.rma(res.rma, 
        xlab=xlab,  slab=sapply(names(esets), .camelCase),
        atransf=exp, at=at, ...)
        return(res.rma)
    }
}

metaCMA.concordance <- function(esets, y="y", risks, rma.method="FE") {
        coefs <- sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~risks[[i]]))$concordance)
        res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        list(res.rma, coefs)
}

metaCMA.hr <- function(esets, y="y", risks, rma.method="FE",inverse=FALSE) {
        coefs <- sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~risks[[i]]))$coefficients[c(1,3)])
        if (inverse) coefs[1,] <- -coefs[1,]
        res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        list(res.rma, coefs)
}

metaCMA.forest <- function(esets, metacma, y="y", mlab="Overall", ...) {
    tmp = names(esets)
    esets = lapply(1:length(esets), function(i) { esets[[i]]$risk =
        metacma$fits[[i]]$risk@lp; esets[[i]]} )   
    names(esets) = tmp    
    .forestplot(esets, y, label="risk", rma.method=metacma$rma.method, mlab=mlab, ...)
}

metaCMA.forest.models <- function(esets, y="y", risks, labeltext=NULL,
summary=TRUE, mlab="Overall",concordance=TRUE,inverse=FALSE,cols=c("darkblue", "seagreen"),...) {
    tmp <- names(esets)
    if (is.null(labeltext)) 
        labeltext <- cbind(c("", sapply(tmp, function(x)
        c(x,rep(NA,length(risks)))),"Overall",rep(NA,length(risks)-1)))
                     #  c("Signature", rep(c("Meta-Analysis","TCGA",NA),
                     #  length(tmp)), c("Meta-Analysis", "TCGA") ))
    if (concordance) {
        rmas <- lapply(risks, function(risk) metaCMA.concordance(esets,y,
            risk))
    } else {
        rmas <- lapply(risks, function(risk) metaCMA.hr(esets,y, risk,
            inverse=inverse))
    }
    r <- do.call(rbind, lapply(1:ncol(rmas[[1]][[2]]),function(i)
    do.call(rbind, c(lapply(rmas, function(rma) rma[[2]][,i]),
        list(c(NA,NA))))))
    r.mean <- c(NA,r[,1], sapply(rmas, function(rma) rma[[1]]$b))
    r.lower <- c(NA,r[,1]-(r[,2]*1.96), sapply(rmas, function(rma) rma[[1]]$ci.lb))
    r.upper <- c(NA,r[,1]+(r[,2]*1.96), sapply(rmas, function(rma)
        rma[[1]]$ci.ub))
    idx <- 1:length(r.mean)
    if (!summary)  idx <- 1:(length(r.mean)-length(risks)-1)

    col=meta.colors(line=c(rep(c(NA, cols),length(tmp)+1))[idx], zero="firebrick",
    box=c(rep(c(NA,cols),length(tmp)+1))[idx])

    forestplot.surv(labeltext=matrix(labeltext[idx,],ncol=ncol(labeltext)), mean=r.mean[idx],
    lower=r.lower[idx],
    upper=r.upper[idx],zero=ifelse(concordance,0.5,0), col=col,
    xlog=concordance==FALSE,
    align=c("l"),  xlab=ifelse(concordance, "Concordance Index", "Hazard Ratio"),
    is.summary=c(rep(FALSE,length(tmp)*(length(risks)+1)+1), rep(TRUE,
    length(risks)))[idx],...)
    rmas
}

metaCMA.forest.probeset <- function(esets, y="y",
probeset, mlab="Overall", rma.method="FE", trans.fun = function(x) x, ...) {
    esets = esets[sapply(esets, function(x) probeset %in% featureNames(x))]
    esets = lapply(esets, function(x) { x[[probeset]] =
    trans.fun(exprs(x)[probeset,]); x })

    .forestplot(esets, y, label=probeset, rma.method=rma.method, mlab=mlab, ...)
}

metaCMA.forest.label <- function(esets, y="y",
label, mlab="Overall", rma.method="FE",...) {
    .forestplot(esets, y, label, rma.method=rma.method, mlab=mlab, ...)
}

metaCMA.common.gene.esets <- function(esets) {
    fns = featureNames(esets[[1]])
    tmp = mapply(function(x) fns <<- intersect(fns,x), lapply(esets, featureNames))
    lapply(esets, function(x) x[fns,])
}

metaCMA.xtable <- function(xtbl, hline.after=getOption("xtable.hline.after", c(-1,0,nrow(xtbl))),
...) {
    print(xtbl,
                  hline.after=NULL,
                  add.to.row=list(pos=lapply(hline.after, function(x) x),
                  command=c("\\toprule\n",
                            rep("\\midrule\n", length(hline.after)-2),
                            "\\bottomrule\n")),...)
}

metaCMA.censor <- function(esets, y="y", censor.at=NULL) {
    .censor <- function(eset) {
        if (!is.null(censor.at)) { 
            eset[[y]][eset[[y]][,1] > censor.at,2] <- 0
            eset[[y]][eset[[y]][,1] > censor.at,1] <- censor.at
        }
        eset
    }
    lapply(esets,.censor)
}

# simple differential expression over all datasets using the limma package
metaCMA.limma <- function(esets, groups, contrasts) {
    require(limma)
    .doLimma <- function(eset) {
        TS <- as.factor(as.character(eset[[groups]]))
        design <- model.matrix(~0+TS)
        colnames(design) <- levels(TS)
        fit <- lmFit(exprs(eset), design = design)
        cont <-
        makeContrasts(contrasts=contrasts,levels=design)
        fit2 <- contrasts.fit(fit, cont)
        fit2 <- eBayes(fit2)
    }    
    lapply(esets, .doLimma)
}

.getCoefsSubset <- function(coefs, idx=1:ncol(coefs[[1]]), probesets=1:nrow(coefs[[1]])) {
    ret <- coefs
    ret[[1]] <- ret[[1]][probesets,idx]
    ret[[2]] <- ret[[2]][probesets,idx]
    ret
}

metaCMA.powerset <- function(n) {
    unlist(lapply(lapply(1:n,function(i) combn(x=1:n, i)), function(x)
        lapply(1:ncol(x), function(i)  x[,i])),recursive=FALSE)
}

metaCMA.allcombinations <- function(i, esets, coefs, eval.fun = metaCMA.eval,
filter.fun=function(eset) return(FALSE), ...) {
    .doPS <- function(ps) {
        idx <- c(i,(1:length(esets))[-i][ps])
        tmp <- metaCMA.train(1, esets[idx],
        coefs=.getCoefsSubset(coefs,idx),filter.fun=filter.fun,...)
        eval.fun(1,esets[idx],object=list(fits=list(tmp)))
     }

     pss <- metaCMA.powerset(length(esets[-i]))
     # use only valid training data (passing our training criteria)
     filtered <- which(sapply(esets[-i], filter.fun))
     pss <- pss[!sapply(pss, function(x) sum(x %in% filtered)>0)]

     ret <- lapply(pss, .doPS)
     n <- lapply(pss, function(ps) sum(sapply(esets[-i][ps], ncol)))
     list(evaluation=ret, n=n)
}

metaCMA.rankproduct <- function(esets, y,
treatment.label=levels(as.factor(esets[[1]][[y]]))[1],
...) {
    require(RankProd)
    esets.c <- .combineEsets(esets, y=y)
    cl <- ifelse(esets.c$y==treatment.label,1,0)
    
    if (sum(cl) == 0 || sum(cl) == length(cl)) stop("Invalid treatment.label.")

    RP.adv.out <- RPadvance(t(esets.c$X),cl,as.numeric(esets.c$batch),...)
}

metaCMA.foldchanges <- function(esets, model, groups,
contrasts) {
    res <- metaCMA.limma(lapply(esets, function(X)
        X[names(model@coefficients),]), groups,
        contrasts)
    resM <- do.call(cbind, lapply(res, function(r)
    topTable(r,number=nrow(esets[[1]]),
        sort.by="none")$logFC))
    resM <- cbind(resM, Weighted.Mean=apply(resM,1, weighted.mean, w=sqrt(sapply(
        esets, ncol))))
    rownames(resM) <- names(model@coefficients)
    resM <- 2^resM
    resM <- apply(resM,c(1,2),function(x) ifelse(x<1, 1/(x*-1),x))
    resM    
}

metaCMA.compare <- function(fe.model, re.model, esets, coefs, ...) {

    fe.all2 <- metaCMA.opt(esets, coefs=coefs, ..., n=nrow(esets[[1]]))
    re.all2 <- metaCMA.opt(esets, coefs=coefs, ..., n=nrow(esets[[1]]),
    rma.method="REML")

    genes <- union(names(fe.model@coefficients), names(re.model@coefficients))
    genes <- genes[genes %in% featureNames(esets[[1]])]

    list(genes = genes, fe = fe.all2, re = re.all2)
}

metaCMA.bootstrap <- function(esets) { 
    lapply(esets, function(X) X[,sample(ncol(X), replace=TRUE)])
}

metaCMA.boxplot <- function(esets, y, probeset, titles) {
    d.f <- do.call(rbind, lapply(1:length(esets), function(i) 
        data.frame(Expression=exprs(esets[[i]])[probeset,], Title=titles[i],
        y=esets[[i]][[y]]))) 
    ggplot(d.f, aes(y,
        Expression))+geom_boxplot()+geom_jitter()+facet_wrap(~Title)
}


