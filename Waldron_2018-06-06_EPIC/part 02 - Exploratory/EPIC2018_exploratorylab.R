## ----setup, echo=FALSE---------------------------------------------------
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(BiocStyle))
suppressPackageStartupMessages(library(curatedMetagenomicData))
# global chunk options
opts_chunk$set(eval=TRUE, cache=TRUE, autodep=TRUE)

## ---- message=FALSE------------------------------------------------------
# BiocInstaller::biocLite("curatedMetagenomicData")
library(curatedMetagenomicData)

## ----eval=FALSE----------------------------------------------------------
## ?combined_metadata
## View(combined_metadata)

## ------------------------------------------------------------------------
table(combined_metadata$antibiotics_current_use)
table(combined_metadata$disease)

## ------------------------------------------------------------------------
dsranking <- combined_metadata %>%
  as.data.frame() %>% 
  group_by(dataset_name) %>%
  summarize(mediandepth = median(number_reads) / 1e6) %>%
  mutate(dsorder = rank(mediandepth)) %>%
  arrange(dsorder)
dsranking

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
combined_metadata %>%
  as.data.frame() %>%
  mutate(ds = factor(combined_metadata$dataset_name, levels=dsranking$dataset_name)) %>%
  ggplot(aes(ds, number_reads / 1e6, fill=body_site)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Dataset", y="Read Depth (millions)")

## ---- message=FALSE------------------------------------------------------
suppressPackageStartupMessages(library(curatedMetagenomicData))
loman.eset = LomanNJ_2013.metaphlan_bugs_list.stool()

## ------------------------------------------------------------------------
loman <- curatedMetagenomicData("LomanNJ_2013.metaphlan_bugs_list.stool", dryrun = FALSE)
loman
loman.eset <- loman[[1]]
loman.eset

## ---- message=FALSE, results='hide'--------------------------------------
oral <- c("BritoIL_2016.metaphlan_bugs_list.oralcavity",
          "Castro-NallarE_2015.metaphlan_bugs_list.oralcavity")
esl <- curatedMetagenomicData(oral, dryrun = FALSE)
esl
esl[[1]]
esl[[2]]

## ---- eval=TRUE----------------------------------------------------------
curatedMetagenomicData("*metaphlan_bugs_list.stool*", dryrun = TRUE)

## ------------------------------------------------------------------------
curatedMetagenomicData("*")

## ------------------------------------------------------------------------
curatedMetagenomicData("*metaphlan*")

## ------------------------------------------------------------------------
curatedMetagenomicData("HMP*")

## ------------------------------------------------------------------------
curatedMetagenomicData("HMP_2012.metaphlan*")

## ------------------------------------------------------------------------
hmp.esetlist = curatedMetagenomicData("HMP_2012.metaphlan*", dryrun=FALSE)
hmp.esetlist

## ------------------------------------------------------------------------
hmp.esetlist = curatedMetagenomicData("HMP_2012.metaphlan*", dryrun=FALSE,
                                     counts=TRUE, bugs.as.phyloseq = TRUE )

## ------------------------------------------------------------------------
eset <- mergeData(esl)
eset

## ------------------------------------------------------------------------
experimentData( loman.eset )

## ------------------------------------------------------------------------
head( pData( loman.eset ) )

## ------------------------------------------------------------------------
exprs( loman.eset )[1:6, 1:5]  #first 6 rows and 5 columns

## ------------------------------------------------------------------------
loman.counts = sweep(exprs( loman.eset ), 2, loman.eset$number_reads / 100, "*")
loman.counts = round(loman.counts)
loman.counts[1:6, 1:5]

## ------------------------------------------------------------------------
loman.eset2 = curatedMetagenomicData("LomanNJ_2013.metaphlan_bugs_list.stool",
                                     counts = TRUE, dryrun = FALSE)[[1]]
all.equal(exprs(loman.eset2), loman.counts)

## ------------------------------------------------------------------------
grep("coli", rownames(loman.eset), value=TRUE)

## ------------------------------------------------------------------------
x = exprs( loman.eset )[grep("s__Escherichia_coli$", rownames( loman.eset)), ]
summary( x )

## ------------------------------------------------------------------------
hist( x, xlab = "Relative Abundance", main="Prevalence of E. Coli",
      breaks="FD")

## ---- warning=FALSE------------------------------------------------------
suppressPackageStartupMessages(library(phyloseq))
loman.pseq = ExpressionSet2phyloseq( loman.eset )

## ------------------------------------------------------------------------
loman.tree <- ExpressionSet2phyloseq( loman.eset, phylogenetictree = TRUE)

## ------------------------------------------------------------------------
wt = UniFrac(loman.tree, weighted=TRUE, normalized=FALSE, 
             parallel=FALSE, fast=TRUE)
plot(hclust(wt), main="Weighted UniFrac distances")

## ------------------------------------------------------------------------
loman.pseq

## ---- warning=FALSE------------------------------------------------------
otu_table( loman.pseq )[1:6, 1:5]

## ---- warning=FALSE------------------------------------------------------
sample_data( loman.pseq )[1:6, 1:5]

## ---- warning=FALSE------------------------------------------------------
head( tax_table( loman.pseq ) )

## ---- warning=FALSE------------------------------------------------------
rank_names( loman.pseq )

## ---- warning=FALSE------------------------------------------------------
subset_taxa( loman.pseq, !is.na(Species))

## ------------------------------------------------------------------------
subset_taxa( loman.pseq, is.na(Class) & !is.na(Phylum))

## ------------------------------------------------------------------------
loman.bd = subset_taxa( loman.pseq, Phylum == "Bacteroidetes")
head( taxa_names( loman.bd ) )

## ---- warning=FALSE------------------------------------------------------
keepotu = genefilter_sample(loman.pseq, filterfun_sample(topp(0.05)), A=5)
summary(keepotu)
subset_taxa(loman.pseq, keepotu)

## ---- warning=FALSE------------------------------------------------------
loman.filt = subset_taxa(loman.pseq, keepotu & !is.na(Strain))
plot_heatmap(loman.filt, method="PCoA", distance="bray")

## ---- warning=FALSE------------------------------------------------------
loman.sp = subset_taxa(loman.pseq, !is.na(Species) & is.na(Strain))
par(mar = c(20, 4, 0, 0) + 0.15) #increase margin size on the bottom
barplot(sort(taxa_sums(loman.sp), TRUE)[1:20] / nsamples(loman.sp),
        ylab = "Total counts", las = 2)

## ---- warning=FALSE------------------------------------------------------
alphas = c("Shannon", "Simpson", "InvSimpson")
plot_richness(loman.sp, "stool_texture", measures = alphas)

## ---- warning=FALSE------------------------------------------------------
pairs( estimate_richness(loman.sp, measures = alphas) )

## ---- warning=FALSE------------------------------------------------------
mydist = distance(loman.sp, method="bray")
myhclust = hclust( mydist )
plot(myhclust, main="Bray-Curtis Dissimilarity", 
     method="ward.D", xlab="Samples", sub = "")

## ---- warning=FALSE------------------------------------------------------
ordinated_taxa = ordinate(loman.sp, method="PCoA", distance="bray")
plot_ordination(loman.sp, ordinated_taxa, color="stool_texture", 
                title = "Bray-Curtis Principal Coordinates Analysis")

## ------------------------------------------------------------------------
plot_scree(ordinated_taxa, title="Screeplot")

