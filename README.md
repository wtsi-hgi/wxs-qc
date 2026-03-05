# WxS-QC - a quality control pipeline for human germline short-variant WGS and WEE cohorts

The **WxS-QC** is a powerful, scalable, and convenient pipeline for 
the QC of human germline short-variant WGS and WES cohorts for population-scale analyses.
The pileine is based on deeply refactored gnomAD v3 and v4 quality control pipelines, 
contains several methods we have developed de novo,
is suitable for both rare-variant discovery and common-variant association studies,
and is aligned with current best practices in WGS/WES germline cohort QC. 

### Pipeline capabilities

We provide:
 
* Well-tested methods from gnomAD v3 and v4.
* Novel QC techniques developed while working with human germline WES and WGS cohorts. In particular:
  - The PCA projection approach that makes PCA-based superpopulation prediction more robust.
  - The “genotype hard filter evaluation” step allows the user to estimate the effect of different filter threshold combinations. 
    Users can choose an optimal combination of variant and genotype-level filters to suit their scientific goal.
* All in one repo, aligned to work together.
* Easily configurable via a YAML-based config file.
* Infrastructure-agnostic solution that can run in any environment supported by the Hail library: a local machine, an on-premises Spark cluster, or a commercial cloud.
* Automatic export of summary tables and interactive graphs.
* A comprehensive, user-friendly tutorial describing the analysis process along with an open test dataset. The complete example of the QC for the public test data is described in the supplementary materials. 
* Excellent scalability due to the usage of the [Hail library](https://hail.is/). The pipeline has been tested on cohorts of 170,000 whole-exome samples and 770 whole-genome samples.
* A reproducible solution that can handle cohorts for different types of experiments in a uniform manner.

### Limitations
* Tested only on short-read data (Illumina, Ultimagenomeics, etc). 
  Technically, there are no restrictions on running the pipeline for long-read data, 
  as long as it provides a standard multi-sample VCF file as input.
* Short-variants: SNVs and small indels. No structural variants support.
* All sample QC approaches are designed for population-scale studies.
  To have robust outlier sample detection, we recommend using a cohort size of >= 150 samples.
* We tested the pipeline extensively only on variants calling by GATK4 suite using haplotype calling and joint calling steps.
  - The pipeline has no direct dependencies of GATK,
    but it uses variant-level and genotype-level metrics calculated by GATK.
    By altering the inputs to the RF model, the variant QC could be adapted for other variant callers 
    (FreeBayes, Strelka2, etc). 
  - The adaptation of the variant QC for modern neural network-based callers, like DeepVariant and DRAGEN, is underway.

## How to QC your data

The high-level description of the QC process is available in a separate document:
[WxS-QC concepts](docs/wxs-qc_concepts.md).

Detailed howto for the QC process is here:
[WcS-QC howto](docs/wxs-qc_howto.md).
For this howto we used an open example dataset, described
[here](docs/wxs-qc_example-dataset.md)

### How to get the latest changes from the `main` branch

When working on an analysis branch, you can retrieve the latest changes from the `main` branch by running:

```bash
make update
```

This will fetch the latest changes from the `main` branch and rebase your current branch onto it.
If there are any unstaged changes in the branch, you will be asked to commit or stash them first.

## Authors

The code is written by members of **Wellcome Sanger HGI group**
(https://www.sanger.ac.uk/group/human-genetics-informatics-hgi/)
based on the **gnomAD QC v3** pipeline by **Broad institute**
(https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v3).


## Developer's howto

For improving the pipeline and developing new functionality,
the [Developer's howto](docs/wxs-qc_development.md) is available.
