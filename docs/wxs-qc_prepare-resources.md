# WxS-QC resources howto

This document describes resource files used by WxS-QC pipeline and how to obtain them.
We expect that you have already set up a work machine / cluster and
cloned the pipeline repository to it, as describen in the
[pipeline howto](wxs-qc_howto.md).

All resources used by WxS-QC pipeline should be stored in the subfolder `resources`
under your analysis folder.

## Downloading WxS-QC resource bundle

The easiest way to get all required resources is to download it from the WSI server:

```bash
wget https://wxs-qc-data.cog.sanger.ac.uk/wxs-qc_resources.tar
tar xvf wxs-qc_resources.tar
```

This archive contains resources from the GATK resource bungle,
the exome set of VCFs for population clustering,
gnomAD population frequencies for exome variations,
and several additional files.
The content of the bundle is described in the details below.

You can copy or symlink the `resources` folder to your data analysis folder.

Please note, that the `resources` folder contains
gnomAD population frequencies and 1000 Genomes data only for exome regions.
To analyze whole-genome data, we suggest using the full-size
gnomAD table For details, see [Using original gnomAD data](#using-original-gnomad-data).
For superpopulation prediction, you can use exome-only 1000 genomes VCFs.


## Resources description

* `1kg_vcf_nosv` -- the latest release of the 1000 Genomes VCFs with removes structural variations.
  For detailed description, see the []()
* `mini_1000G` -- one chromosome from the `1kg_vcf_nosv` with reduced number of samples.
  Used for test purposes. Usually, you don't need to use it.
* `igsr_samples.tsv` -- known super populations for 1000 genomes dataset.
  Available here: https://www.internationalgenome.org/data-portal/sample (press the 'Download the list' button)
* `long_ld_regions.hg38.bed` -- BED file containing long-range linkage disequilibrium regions for the genome version hg38
  The regions were obtained from the file `high-LD-regions-hg38-GRCh38.bed` in **plinkQC** github repo:
  (https://github.com/cran/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed).
  These coordinates are results of `liftOver` transferring original coordinates from the genome version hg36 to hg38.
  Original coordinates are provided in supplementary files of the article
  **Anderson, Carl A., et al. "Data quality control in genetic case-control association studies."
  Nature protocols 5.9 (2010): 1564-1573. DOI: 10.1038/nprot.2010.116**
* `HG001_GRCh38_benchmark.interval.illumina.vcf.gz` with `tabix` index -- High-confidence variations for GIAB HG001 sample.
  Available here: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
* `HG001_GRCh38_benchmark.all.interval.illumina.vep.info.txt` - VEP annotations for GIAB HG001 sample.
* `gnomad.exomes.r4.1.freq_only.ht` - reduced version of **gnomAD** data containing only global population frequencies
* `1000G_phase1.snps.high_confidence.hg38`, `1000G_omni2.5.hg38`,
  `hapmap_3.3.hg38`, `Mills_and_1000G_gold_standard.indels.hg38` -
   the set of high-confident variations from the corresponding projects.
   The corresponding VCFs are available from the GATK resource bundle:
   (https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).
   - Google Cloud: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
   - FTP: http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg38/


## Building 1000 Genomes VCFs

The SampleQC stage of the WxS-QC pipeline uses data from
[1000Genomes VCF dataset](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/),
to run clustering and identify superpopulations.
This section describes how this data was made.

### Prepare a work folder

Make the directory for 1000 genomes data

```bash
export KG_DIR=/path/to/your/dataset/folder/resources/1kg_vcf
mkdir -p $KG_DIR
cd $KG_DIR
```

### Download data

The `1kg_download` script downloads VCF files from the 1000 Genomes Project:

```bash
/wes/qc/root/folder/scripts/1kg_download
```

Make a symlink the VCF for X-chromosome, because for some reason it has a slightly different name:

```bash
ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
```


### Remove structural variants
The downloaded 1000 genomes VCFs contain structural variation.
We remove it and keep only SNVs and short indels, that are used for PCA and clustering.

The `1kg_remove_sv` script removes structural variants from VCF files.
You need to run the script from the 1000 genomes data folder with the chromosome name as argument.

```bash
# To process all chromosomes sequentially
for chr in {1..22} X; do /wes/qc/root/folder/scripts/1kg_remove_sv $chr; done
```

If you work on a cluster with the installed IBM LSF job scheduling system,
you can use it to submit your jobs.
The script automatically uses LSF job array number for chromosome names

```bash
# Chromosome X
bsub -n 2 /wes/qc/root/folder/scripts/1kg_remove_sv X

# Pick all numeric chromosomes from the job array ID
bsub -n 2 -J [1-22] /wes/qc/root/folder/scripts/1kg_remove_sv
```

### Remove original files

After checking that all files were processed correctly,
you can remove the original 1000 genomes files:

```bash
rm 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

### (Optional) Generate exome-only VCFs

The manual above prepares the full-genome dataset, containing variations across the whole human genome,
If you work only with exome data, you can significantly reduce the data size,
keeping only exome regions.

First, you need the BED/GFF/GTF file containing coordinates of exome regions.

If you have the BED file for your exome enrichment kit,
we highly recommend you to use it.

If you don't have any, you can download GTF from _Gencode_ and keep only exon regions from it.

First, download all reference data and create genome size file from reference genome.
You need the `samtools` package to generate index containing genome sizes.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

gunzip GRCh38.primary_assembly.genome.fa.gz
samtools faidx GRCh38.primary_assembly.genome.fa
```

Now extract only exons, add +50 bp to capture flanking regions and merge overlapping regions.
For this work you need the `bedtools` package to manipulate genome regions.

```bash
gunzip -c gencode.v47.annotation.gtf.gz | awk '{FS="\t";OFS="\t"} $3=="exon"' | bedtools slop -g GRCh38.primary_assembly.genome.fa.fai -b 50 | bedtools sort -g GRCh38.primary_assembly.genome.fa.fai | bedtools merge > GRCh38.gencode.v47.exome.bed
```

Finally, extract only exome regions from 1000 genomes VCFs:

```bash
for chr in {1..22} X; do bcftools view -R GRCh38.gencode.v47.exome.bed 1kGP_chr${chr}.nosv.vcf.gz -Oz > 1kGP_chr${chr}.nosv.exome.vcf.gz; done
```

Now you can copy all VCFs with the `nosv.exome.vcf.gz` suffix to a new folder,
and use it in step 0.1 of the resource preparation stage (stage 0).

## Using original gnomAD data

The Variant QC part of the pipeline uses population frequencies from the
[gnomAD project](https://gnomad.broadinstitute.org/)
to find _de novo_ variations.
Technically, for this step you can use the original gnomAD exome/genome data.
However, the full-size gnomAD dataset is very big, so we prepared
a reduced version, containing only global population frequencies for exome regions.

If you want to use your own gnomAD data (for example, for whole-genome data),
  you need to manually download it from https://gnomad.broadinstitute.org/downloads
  (use the _Sites Hail Table_ version),
  place the path to the table in the config file section `prepare_gnomad_ht -> input_gnomad_htfile`,
  and run the script to make a reduced version:
  ```shell
  spark-submit 0-resource_preparation/3-prepare-gnomad-table.py
  ```
