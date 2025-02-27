# RSV-Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/build_using-nf--core-1a9655)](https://nf-co.re/)

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Follow on X](http://img.shields.io/badge/%40LabCFlores-1DA1F2?labelColor=000000&logo=X)](https://x.com/LabCFlores)

A public repository of **Respiratory Syncytial Virus genomic surveillance** bioinformatic pipeline maintained by ITER.

## Introduction

**RSVPipeline** automates the processing of Illumina and Oxford Nanopore Technologies (ONT) sequencing datasets for viral genome analysis. It integrates multiple tools for quality control, taxonomic classification, reference selection, alignment, consensus generation, and variant calling.

The **RSVPipeline** is built using [Nextflow](https://www.nextflow.io/), following [nf-core](https://nf-co.re) guidelines and templates.

## Illumina pipeline summary

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/Illumina_pipeline_GitHub.jpg">
    <img alt="Illumina pipeline" src="docs/images/Illumina_pipeline_GitHub.jpg">
  </picture>
</h1>

1. Quality control of raw reads with FastQC
2. Taxonomic classification of raw reads using Kraken2 with PlusPF database
3. Adapter trimming with fastp
4. Quality control of trimmed reads with FastQC
5. Remove host reads using Kraken2 with humanDB database
6. Multi-reference alignment with bbmap to select reference
  6.1. Run alternative assembly strategy if previous MSA step failed
7. Align reads to the correct reference strain
8. Trim adapters with iVar before create consensus sequences
9. Coverage analysis with Mosdepth
10. Create consensus sequence with iVar
11. Variant-calling with iVar

## ONT pipeline summary

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/ONT_pipeline_GitHub.jpg">
    <img alt="ONT pipeline" src="docs/images/ONT_pipeline_GitHub.jpg">
  </picture>
</h1>

1. Quality control of raw reads with NanoPlot
2. Taxonomic classification of raw reads using Kraken2 with PlusPF database
3. Remove host reads using Kraken2 with humanDB database
4. Quality control of dehosted reads with NanoPlot
5. Select reference using IRMA
6. Align reads to the correct reference strain
7. Trim adapters with align_trim (from ARTIC) before create consensus sequences
8. Coverage analysis with Mosdepth
9. Create consensus sequence with Medaka
10. Variant-calling with Medaka

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.





## How to cite this work

This work has not been publised yet.

Please cite this repository as: "RSVPipeline (accessed on YYYY-MM-DD)". And do not forget to cite the paper when it becomes available.

## Funding

...
