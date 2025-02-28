# nf-rsvpipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/build_using-nf--core-1a9655)](https://nf-co.re/)

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Follow on X](http://img.shields.io/badge/%40LabCFlores-1DA1F2?labelColor=000000&logo=X)](https://x.com/LabCFlores)

A public repository of **Respiratory Syncytial Virus genomic surveillance** bioinformatic pipeline maintained by ITER.

## Introduction

**nf-rsvpipeline** automates the processing of Illumina and Oxford Nanopore Technologies (ONT) **amplicon-based** sequencing datasets for viral genome analysis. It integrates multiple tools for quality control, taxonomic classification, reference selection, alignment, consensus generation, and variant calling.

The **nf-rsvpipeline** is built using [Nextflow](https://www.nextflow.io/), following [nf-core](https://nf-co.re) guidelines and templates.

### Respiratory Syncytial Virus genomic surveillance in the Canary Islands

This pipeline is designed for the analysis of RSV whole genomes using short- and long-read sequencing technologies. It is developed as part of the research efforts documented in the following repository:

[`RSV Repository (genomicsITER/RSV)`](https://github.com/genomicsITER/RSV)

The RSV repository contains genomic data, analysis scripts, and additional resources related to the study of Respiratory Syncytial Virus (RSV) cases detected in the Canary Islands between 2022 and 2024. The pipeline in this repository facilitates the automated processing and variant analysis of RSV sequencing datasets.

## Illumina pipeline summary

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/Illumina_pipeline_GitHub.jpg">
    <img alt="Illumina pipeline" src="docs/images/Illumina_pipeline_GitHub.jpg">
  </picture>
</h1>

1. Quality control of raw reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Taxonomic classification of raw reads ([`Kraken2`](https://github.com/DerrickWood/kraken2))
   * Use PlusPF database ([`Kraken2 databases`](https://benlangmead.github.io/aws-indexes/k2))
3. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp))
4. Quality control of trimmed reads ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
5. Remove host reads ([`Kraken2`](https://github.com/DerrickWood/kraken2))
   * Use HumanDB ([`Kraken2 databases`](https://benlangmead.github.io/aws-indexes/k2))
6. Multi-reference alignment to select reference ([`BBMap`](https://sourceforge.net/projects/bbmap/))
   1. Run alternative assembly strategy if previous MSA step failed ([`SPAdes`](https://github.com/ablab/spades))
7. Align reads to the correct reference strain ([`BWA`](https://github.com/lh3/bwa/))
8. Trim adapters before create consensus sequences ([`iVar`](https://github.com/andersen-lab/ivar))
9. Coverage analysis ([`MosDepth`](https://github.com/brentp/mosdepth))
10. Create consensus sequence ([`iVar`](https://github.com/andersen-lab/ivar))
11. Variant-calling ([`iVar`](https://github.com/andersen-lab/ivar))

## ONT pipeline summary

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/ONT_pipeline_GitHub.jpg">
    <img alt="ONT pipeline" src="docs/images/ONT_pipeline_GitHub.jpg">
  </picture>
</h1>

1. Quality control of raw reads ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
2. Taxonomic classification of raw reads ([`Kraken2`](https://github.com/DerrickWood/kraken2))
   * Use PlusPF database ([`Kraken2 databases`](https://benlangmead.github.io/aws-indexes/k2))
3. Remove host reads ([`Kraken2`](https://github.com/DerrickWood/kraken2))
   * Use HumanDB ([`Kraken2 databases`](https://benlangmead.github.io/aws-indexes/k2))
4. Quality control of dehosted reads ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
5. Select reference to downstream analysis ([`IRMA`](https://wonder.cdc.gov/amd/flu/irma/index.html))
6. Align reads to the correct reference strain ([`Minimap2`](https://github.com/lh3/minimap2))
7. Trim adapters before create consensus sequences ([`ARTIC`](https://github.com/artic-network/fieldbioinformatics))
8. Coverage analysis ([`MosDepth`](https://github.com/brentp/mosdepth))
9. Create consensus sequence ([`Medaka`](https://github.com/nanoporetech/medaka))
10. Variant-calling ([`Medaka`](https://github.com/nanoporetech/medaka))

## Quick start

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=24.04.2`) before running the pipeline.

2. Make sure to download and set up the **Kraken2 PlusPF database** from this [`repository`](https://benlangmead.github.io/aws-indexes/k2). Ensure that the database path is correctly set in the `nextflow.config` file.

3. Clone the repository:

   ```bash
   git clone https://github.com/genomicsITER/nf-rsvpipeline
   cd nf-rsvpipeline
   ```

4. Input Data Requirements: To run the pipeline, the user must provide sequencing data in the following format:

   * **For Illumina data**: A folder containing paired-end FASTQ files (`_R1.fastq.gz` and `_R2.fastq.gz` for each sample).
   * **For ONT data**: A folder containing single-end FASTQ files (`.fastq.gz`).

5. Run the pipeline:

   For Illumina **paired-end** datasets run the pipeline as follows:

   ```bash
   nextflow run main.nf \
     --indir /path/to/fastq_dir/illumina \
     --outdir results_illumina \
     --platform "illumina" \
     -profile <docker/singularity/conda> \
     -resume
   ```

   For ONT datasets run the pipeline as follows:

   ```bash
   nextflow run main.nf \
     --indir /path/to/fastq_dir/nanopore \
     --outdir results_nanopore \
     --platform "nanopore" \
     -profile <docker/singularity/conda> \
     -resume
   ```

   Additionaly, if you want to run [`pycoQC`](https://github.com/a-slide/pycoQC), you need to add the `--sequencing_summary` parameter:

   ```bash
   nextflow run main.nf \
     --indir /path/to/nanopore/fastq/files \
     --outdir results_nanopore \
     --platform "nanopore" \
     --sequencing_summary sequencing_summary.txt \
     -profile <docker/singularity/conda> \
     -resume
   ```

For further assistance, feel free to open an [issue](https://github.com/genomicsITER/nf-rsvpipeline/issues) in this repository.

## Pipeline output

Overview of the different results produced by the pipeline and how to interpret them.

<details>

<summary>**Summary output of Illumina pipeline**</summary>

### Preprocessing

* `01_FastQC`
  * `*_{1,2}_fastqc.html`: FastQC report containing quality metrics.
  * `*_{1,2}_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
* `02_`
  * `*.kraken2.classifiedreads.txt`: .
  * `*.kraken2.report.txt`: .
* `03_Adapter_trimming`
  * `*`: .
  * `*.fastp.log`: .
  * `*.fastp.json`: .
  * `*.fastp.html`: .
  * `*_{1,2}.fastp.fastq.gz`: .
* `04_FastQC_trimmed`
  * `*_{1,2}_fastqc.html`: FastQC report containing quality metrics of trimmed resds.
  * `*_{1,2}_fastqc.zip`: Zip archive containing the FastQC-trimmed report, tab-delimited data file and plot images.
* `05_Remove_host_reads`
  * `*.classified_{1,2}.fastq.gz`: .
  * `*.unclassified_{1,2}.fastq.gz`: .
  * `*.kraken2.classifiedreads.txt`: .
  * `*.kraken2.report.txt`: .
* `06_Reference_selection`
  * `*.MSA.bam`: .
  * `*.MSA.bbmap.log`: .
  * `*.MSA.covstats.tsv`: .
  * `*_refs.tsv`: .
* `07_Align`
  * `*.RSV-{A,B}.bam`: .
  * `*.RSV-{A,B}.bam.bai`: .
  * `*.RSV-{A,B}.flagstat`: .
  * `*.RSV-{A,B}.idxstats`: .
  * `*.RSV-{A,B}.mapped.bam`: .
  * `*.RSV-{A,B}.mapped.bam.bai`: .
* `08_Trim_adapters`
  * `*.RSV-{A,B}.mapped.trimmed.bam`: .
  * `*.RSV-{A,B}.mapped.trimmed.ivar.log`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.bam`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.bam.bai`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.flagstat`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.idxstats`: .
* `09_Coverage_analysis`
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.mosdepth.global.dist.txt`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.mosdepth.region.dist.txt`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.mosdepth.summary.txt`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.per-base.bed.gz`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.per-base.bed.gz.csi`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.regions.bed.gz`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.all_amplicons.regions.bed.gz.csi`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.mosdepth.global.dist.txt`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.mosdepth.summary.txt`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.per-base.bed.gz`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.per-base.bed.gz.csi`: .
  * `*.RSV-{A,B}.mapped.trimmed.sorted.tsv`: .
* `10_Calling`
  * `*.RSV-{A,B}.ivar_calling.mpileup`: .
  * `*.RSV-{A,B}.ivar_calling.tsv`: .
  * `*.RSV-{A,B}.ivar_calling.vcf`: .
* `11_Consensus`
  * `*.RSV-{A,B}.ivar_consensus.fa`: .
  * `*.RSV-{A,B}.ivar_consensus.mpileup`: .
  * `*.RSV-{A,B}.ivar_consensus.qual.txt`: .

You can add text within a collapsed section.

You can add an image or a code block, too.

```ruby
   puts "Hello World"
```

</details>

## How to cite this work

This work has not been publised yet.

Please cite this repository as: "nf-rsvpipeline (accessed on YYYY-MM-DD)". And do not forget to cite the paper when it becomes available.

## Funding

...
