# WxS-QC concepts

The aim of the WxS-QC pipeline is to produce a clean dataset removing problem samples  
and excluding variants and genotypes that are likely to be false positives. 

Each dataset is different,
and the QC requires careful tailoring and evaluation at each step
to ensure that it is appropriate for the dataset being worked on.

## Overview

WxS-QC consists of four main steps: data loading, sample QC, variant QC, and genotype QC. 

It is written using [Hail](https://hail.is/) 
and much of it is based on methods used by gnomAD:
(https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/,   
https://gnomad.broadinstitute.org/news/2019-10-gnomad-v3-0/,  
https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/)

### Data loading and external validation

The first block converts raw VCFs to Hail data structures and validates them against the provided external metadata:
* Checks the threshold for VerifyBamID FreeMix score
* If the microarray genotype data are available, the pipeline checks consistency between sequencing calls 
  and microarray data and flags potentially problematic samples.

### Sample QC

The aim of sample QC is to remove poor-quality samples.
This block generally follows the approach described in 
(Sealock et al. _"Tutorial: guidelines for quality filtering of whole-exome and whole-genome sequencing data 
for population-scale association analyses"_. Nat Protoc 2025:1–11.):
make a high-quality variation subset, impute sex,
check sex consistency, and apply the stratified threshold to filter out outlier samples.
We implemented three different methods for outlier detection:
  * `pop` - Stratified QC with PCA-assigned superpopulations with PCA projection (from gnomAD v3)
  * `nn` - Nearest neighbors (from gnomAD v4)
  * `lr` - Linear regression (from gnomAD v4)

For details, eee the [Sample QC stage documentation](wxs-qc_howto.md#stage-2-sample-qc).

## Variant QC

The VariantQC step uses the approach used for gnomAD v2 (Karczewski et al. 2020) and the first stages of gnomAD v3.

We use a set of open resources to label variations as likely True-Positives (TP) and likely False-Positives (FP). 
Then we train a random forest (RF) model to predict TP and FP variants based on variant-level statistics, 
use RF model score to group variants into several bins based on their reliability.

Historically, this step was designed for “classic” variant callers, which require recalibrating the call correctness. 
In the development version, we’re updating this part to support DeepVariant 
and other neural network-based callers, which provide a very brief set of variant-level statistics.

We have tested this extensively only on [GATK v4](https://gatk.broadinstitute.org/hc/en-us) HaplotypeCaller. 
However, by altering the inputs to the RF model, 
the variant QC could be adapted for other variant callers (FreeBayes, Strelka2, etc).

## Genotype QC

The final module, hard filter evaluation and genotype QC, is a new functionality developed at WSI.
The detailed description of this step is available in the pipeline documentation and in the paper (Koko et al. 2024).
At this step,
we test
how different combinations of variant-level and genotype-level filters affect the genotypes and variations surviving,
plot and review the results, and choose optimal combinations.

