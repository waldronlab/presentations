## ----setup, cache=FALSE, echo=FALSE--------------------------------------
library(knitr)
# global chunk options
opts_chunk$set(cache=TRUE, autodep=TRUE)

## ---- echo=FALSE---------------------------------------------------------
par(mar=c(4, 4, 0, 0))
plot(x=0:10, y=dpois(0:10, lambda=1), 
     type="b", lwd=2,
     xlab="Counts (k)", ylab="Probability density")
lines(x=0:10, y=dpois(0:10, lambda=2), 
      type="b", lwd=2, lty=2, pch=2)
lines(x=0:10, dpois(0:10, lambda=4), 
      type="b", lwd=2, lty=3, pch=3)
legend("topright", lwd=2, lty=1:3, pch=1:3,
       legend=c(expression(paste(lambda, "=1")),
                expression(paste(lambda, "=2")),
                expression(paste(lambda, "=4"))))

## ---- echo=FALSE---------------------------------------------------------
plot(x=0:40, y=dnbinom(0:40, size=10, prob=0.5), 
     type="b", lwd=2, ylim=c(0, 0.2),
     xlab="Counts (k)", ylab="Probability density")
lines(x=0:40, y=dnbinom(0:40, size=20, prob=0.5), 
      type="b", lwd=2, lty=2, pch=2)
lines(x=0:40, y=dnbinom(0:40, size=10, prob=0.3),
      type="b", lwd=2, lty=3, pch=3)
legend("topright", lwd=2, lty=1:3, pch=1:3,
       legend=c("n=10, p=0.5", "n=20, p=0.5", "n=10, p=0.3"))

## ---- echo=FALSE---------------------------------------------------------
plot(x=0:40, y=dnbinom(0:40, size=10, prob=0.5), 
     type="b", lwd=2, ylim=c(0, 0.15),
     xlab="Counts (k)", ylab="Probability density")
lines(x=0:40, y=dnbinom(0:40, size=20, prob=0.5), 
      type="b", lwd=2, lty=2, pch=2)
lines(x=0:40, y=dnbinom(0:40, size=10, prob=0.3),
      type="b", lwd=2, lty=3, pch=3)
lines(x=0:40, y=dpois(0:40, lambda=9), col="red")
lines(x=0:40, y=dpois(0:40, lambda=20), col="red")
legend("topright", lwd=c(2,2,2,1), lty=c(1:3,1), pch=c(1:3,-1), col=c(rep("black", 3), "red"),
       legend=c("n=10, p=0.5", "n=20, p=0.5", "n=10, p=0.3", "Poisson"))

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library(gamlss)
##par(cex=2)  #increase size of type and axes
plot(x=0:10, y=dpois(0:10, lambda=2), 
     type="b", lwd=2, ylim=c(0, 0.5),
     xlab="Counts (k)", ylab="Probability density")
lines(x=0:10, y=dZIP(0:10, mu=2, sigma=0.2),
      type="b", lwd=2, lty=2, pch=2)
lines(x=0:10, y=dZIP(0:10, mu=2, sigma=0.4),
      type="b", lwd=2, lty=3, pch=3)
legend("topright", lwd=2, lty=1:3, pch=1:3,
       legend=c(expression(paste(lambda, "=2")),
                expression(paste("ZIP: ", lambda, "=2, ", "p=0.2")),
                expression(paste("ZIP: ", lambda, "=2, ", "p=0.4"))))

## ------------------------------------------------------------------------
indat = read.delim("data/Candela_Africa_stool.txt")
inmetadat = read.delim("data/Candela_Africa_metadat.txt")
library(phyloseq)
source("https://raw.githubusercontent.com/waldronlab/EPIC2016/master/metaphlanToPhyloseq.R")
Candela = metaphlanToPhyloseq(indat, inmetadat, simplenames = TRUE, 
                              roundtointeger = FALSE)

## ---- eval=FALSE---------------------------------------------------------
## summary(otu_table(Candela))
## summary(sample_data(Candela))

## ------------------------------------------------------------------------
Candela
Candela = prune_taxa(taxa_sums(Candela) > 0, Candela)
Candela

## ------------------------------------------------------------------------
rank_names(Candela)
subset_taxa(Candela, !is.na(Strain))
(Candela.sp_strain = subset_taxa(Candela, !is.na(Species)))
taxonomy.level = apply(tax_table(Candela), 1, function(x) sum(!is.na(x)))
Candela.phy = prune_taxa(taxonomy.level==2, Candela)
taxa_names(Candela.phy)

## ------------------------------------------------------------------------
f1<- filterfun_sample(topp(0.1))
pru <- genefilter_sample(Candela, f1, A=2)
summary(pru)
subset_taxa(Candela, pru)

## ------------------------------------------------------------------------
plot_heatmap(Candela.sp_strain, method="PCoA", distance="bray")

## ------------------------------------------------------------------------
par(mar = c(18, 4, 0, 0) + 0.1) # make more room on bottom margin
barplot(sort(taxa_sums(Candela.sp_strain), TRUE)[1:30]/nsamples(Candela.sp_strain), las=2)

## ---- echo=FALSE, fig.height=3.5-----------------------------------------
rafalib::mypar()
plot(c(0,1,1),c(0,0,1),pch=16,cex=2,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(-0.25,1.25),ylim=c(-0.25,1.25))
lines(c(0,1,1,0),c(0,0,1,0))
text(0,.2,expression(paste('(A'[x]*',A'[y]*')')),cex=1.5)
text(1,1.2,expression(paste('(B'[x]*',B'[y]*')')),cex=1.5)
text(-0.1,0,"A",cex=2)
text(1.1,1,"B",cex=2)

## ---- warning=FALSE------------------------------------------------------
Candela.int = Candela
otu_table(Candela.int) = round(otu_table(Candela)*1e6)
alpha_meas = c("Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(Candela.int, "gender", "camp", measures=alpha_meas))

## ------------------------------------------------------------------------
alphas = estimate_richness(Candela.int, measures=alpha_meas)
pairs(alphas)

## ------------------------------------------------------------------------
plot(hclust(phyloseq::distance(Candela, method="bray")), 
     main="Bray-Curtis Dissimilarity", xlab="", sub = "")

## ------------------------------------------------------------------------
ord = ordinate(Candela, method="PCoA", distance="bray")
plot_ordination(Candela, ord, color="camp", shape="camp") + 
  ggplot2::ggtitle("Bray-Curtis Principal Coordinates Analysis")

## ------------------------------------------------------------------------
plot_scree(ord) + ggplot2::ggtitle("Screeplot")

## ---- eval=FALSE---------------------------------------------------------
## Candela.pruned = subset_samples(Candela, country %in% c("italy", "tanzania"))
## Candela.pruned = prune_samples(sample_sums(Candela.pruned) > 10, Candela.pruned)
## Candela.pruned

## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library("DESeq2"))

## ------------------------------------------------------------------------
dds.data = phyloseq_to_deseq2(Candela, ~country)

## ------------------------------------------------------------------------
dds = DESeq(dds.data)
res = results(dds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Candela)[rownames(sigtab), ], "matrix"))
head(sigtab)

## ------------------------------------------------------------------------
library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Family))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Family order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
ggplot(sigtabgen, aes(y=Family, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## ------------------------------------------------------------------------
table(sample_data(Candela)$country, sample_data(Candela)$gender)

## ------------------------------------------------------------------------
dds.data2 = phyloseq_to_deseq2(Candela, ~country + gender)
dds2 = DESeq(dds.data)
resultsNames(dds2)

## ------------------------------------------------------------------------
## italy = numerator, tanzania = denominator
res2 = results(dds, contrast=c("country", "italy", "tanzania"))
res2 = res2[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab2 = res2[which(res2$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(Candela)[rownames(sigtab2), ], "matrix"))
head(sigtab2)

## ------------------------------------------------------------------------
plotMA(res, main="Difference vs. Average")
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
plotCounts(dds2, gene="p__Actinobacteria", intgroup="country")
plotCounts(dds2, gene="p__Actinobacteria", intgroup="gender")

## ------------------------------------------------------------------------
select <- rownames(sigtab2)
nt <- normTransform(dds2) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select, ]
df <- as.data.frame(colData(dds2)[,c("country", "gender")])

## ------------------------------------------------------------------------
pheatmap::pheatmap(log2.norm.counts, annotation_col=df, main="log2(counts + 1)")

## ------------------------------------------------------------------------
df = data.frame(country=sample_data(Candela)$country,
                Shannon=alphas$Shannon)
df = cbind(df, ord$vectors[, 1:5])

## ------------------------------------------------------------------------
par(mfrow=c(3,2))
for (i in 2:7){
  boxplot(df[, i] ~ df$country, main=colnames(df)[i])
}

## ------------------------------------------------------------------------
fit = glm(country ~ ., data=df, family=binomial("logit"))

## ------------------------------------------------------------------------
res <- sapply(2:ncol(df), function(i){
  fit = glm(df$country ~ df[, i], family=binomial("logit"))
  summary(fit)$coefficients[2, ]
})
colnames(res) = colnames(df)[2:ncol(df)]
res
write.csv(res, file="univariate_shannonPCoA.csv")

