# Public Data Resources and Bioconductor

# Instructors

- Levi Waldron, City University of New York, New York, NY, USA
- Benjamin Haibe-Kain, Princess Margaret Cancer Center, Toronto, Canada
- Sean Davis, Center for Cancer Research, National Cancer Institute, National Institutes of Health, Bethesda, MD, USA

# Workshop Description 

The goal of this workshop is to introduce Bioconductor packages for finding,
accessing, and using large-scale public data resources including the 
Gene Expression Omnibus [GEO](https://www.ncbi.nlm.nih.gov/geo), Sequence
Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra), the Genomic Data
Commons [GDC](https://portal.gdc.cancer.gov/), and Bioconductor-hosted 
curated data resources for metagenomics, pharmacogenomics, and The Cancer 
Genome Atlas.

## Pre-requisites

* Basic knowledge of R syntax
* Familiarity with the ExpressionSet and SummarizedExperiment classes
* Basic familiarity with 'omics technologies such as microarray and NGS sequencing

Interested students can prepare by reviewing vignettes of the following packages to gain background on aspects of interest to them:

* [GEOquery](http://bioconductor.org/packages/GEOquery/): Access to the NCBI Gene Expression Omnibus (GEO), a public repository of gene expression (primarily microarray) data.
* [GenomicDataCommons](http://bioconductor.org/packages/GenomicDataCommons/): Access to the NIH / NCI Genomic Data Commons RESTful service.
* [SRAdb](http://bioconductor.org/packages/SRAdb/): A compilation of metadata from the NCBI Sequence Read Archive, the largest public repository of sequencing data from the next generation of sequencing platforms, and tools
* [curatedTCGAData](http://bioconductor.org/packages/curatedTCGAData/): Curated data from The Cancer Genome Atlas (TCGA) as MultiAssayExperiment Objects
* [curatedMetagenomicData](http://bioconductor.org/packages/curatedMetagenomicData/): Curated metagenomic data of the human microbiome
* [PharmacoGx](https://bioconductor.org/packages/PharmacoGx/): Analysis of large-scale pharmacogenomic data

## workshop Participation 

Describe how students will be expected to participate in the workshop.

## R/Bioconductor packages used

List any R/Bioconductor packages that will be explicitly covered.

## Time outline

An example for a 45m workshop:

| Activity                            | Time    |
|-------------------------------------|---------|
| Packages | 15m |
| Package Development | 15m |
| Contributing to Bioconductor | 5m |
| Best Practices | 10m |

# workshop goals and objectives

List "big picture" student-centered workshop goals and learning
objectives. Learning goals and objectives are related, but not the
same thing. These goals and objectives will help some people to decide 
whether to attend the conference for training purposes, so please make 
these as precise and accurate as possible.

*Learning goals* are high-level descriptions of what
participants will learn and be able to do after the workshop is
over. *Learning objectives*, on the other hand, describe in very
specific and measurable terms specific skills or knowledge
attained. The [Bloom's Taxonomy](#bloom) may be a useful framework 
for defining and describing your goals and objectives, although there
are others.

## Learning goals

Some examples:

* describe how to...
* identify methods for...
* understand the difference between...

## Learning objectives

* analyze xyz data to produce...
* create xyz plots
* evaluate xyz data for artifacts



### A note about learning goals and objectives (#bloom)

While not a new or modern system for thinking about learning,
[Bloom's taxonomy][1] is one useful framework for understanding the
cognitive processes involved in learning. From lowest to highest
cognitive requirements:

1. Knowledge: Learners must be able to recall or remember the
   information.
2. Comprehension: Learners must be able to understand the information.
3. Application: Learners must be able to use the information they have
   learned at the same or different contexts.
4. Analysis: Learners must be able to analyze the information, by
   identifying its different components.
5. Synthesis: Learners must be able to create something new using
   different chunks of the information they have already mastered. 
6. Evaluation: Learners must be able to present opinions, justify
   decisions, and make judgments about the information presented,
   based on previously acquired knowledge.
   
To use Bloom's taxonomy, consider the following sets of verbs and
descriptions for learning objectives:

1. Remember: Memorize, show, pick, spell, list, quote, recall, repeat,
   catalogue, cite, state, relate, record, name.
2. Understand: Explain, restate, alter, outline, discuss, expand,
   identify, locate, report, express, recognize, discuss, qualify,
   covert, review, infer.
3. Apply: Translate, interpret, explain, practice, illustrate,
   operate, demonstrate, dramatize, sketch, put into action, complete,
   model, utilize, experiment, schedule, use.
4. Analyze: Distinguish, differentiate, separate, take apart,
   appraise, calculate, criticize, compare, contrast, examine, test,
   relate, search, classify, experiment.
5. Evaluate: Decide, appraise, revise, score, recommend, select,
   measure, argue, value, estimate, choose, discuss, rate, assess,
   think.
6. Create: Compose, plan, propose, produce, predict, design, assemble,
   prepare, formulate, organize, manage, construct, generate, imagine,
   set-up.

[1]: https://cft.vanderbilt.edu/guides-sub-pages/blooms-taxonomy/ "Bloom's Taxonomy"
