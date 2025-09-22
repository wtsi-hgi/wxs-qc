# The WxS-QC example dataset
The public example dataset we use to test the WxS-QC pipeline is prepared from the open 1000 genomes data,
and is freely available for download from https://wxs-qc-data.cog.sanger.ac.uk/wxs-qc_public_dataset_v3.tar.
This dataset contains the combined VCF file `joint_germline.vcf.gz` and the set of metadata files:
* `verify_bam_id_result.selfSM.tsv` - VerifyBamID results
* `self_reported_sex.tsv` - self-reported sex for all 1000 genomes samples.
* `trios.fam`, `trios-withparents.fam` – pedigree information. The pipeline now uses the version without parents.
  Switching to the full version with parents is in progress.
* `all_consequences_with_gene_and_csq.tsv.bgz`, `all_consequences.tsv.bgz`, `synonymous_variants.tsv.bgz` -
  VEP annotation in different formats. For now, the pipeline uses it all. The work on unifying the VEP input is in progress.

In the dataset-build subfolder, you can find all the scripts we used to prepare the data.
The scripts are numbered and contain all the details to trace data preparation processes and ensure reproducibility.

Briefly, the preparation process is the following:
* Download the list of high coverage CRAMS from
  [1000 genomes FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/),
  released 2015-12-22. The list of files contains a pre-selected set of samples to cover all super populations
  and mimic outliers in real experimental datasets.
* Run [NF-Core Sarek pipeline](https://github.com/nf-core/sarek) with GATK4 HaplotypeCaller
  to call per-sample variations and make joint-calling. The resulting VCF is the initial data for the QC pipeline.
* Run VerifyBamID on CRAMs and collect FreeMix scores.
* Download the file with self-reported sex information and change the header to match the pipeline requirements.
* Run VEP on the VCF to obtain consequences. To speed up calculations, we normalize the VCF and split it into several chunks.
