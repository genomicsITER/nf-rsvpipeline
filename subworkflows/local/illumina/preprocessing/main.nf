//
// Subworkflow with preprocessing functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC as FASTQC_RAW               } from '../../../../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_PLUSPF  } from '../../../../modules/nf-core/kraken2/kraken2/main'
include { FASTP                              } from '../../../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_TRIM              } from '../../../../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HUMANDB } from '../../../../modules/nf-core/kraken2/kraken2/main'
include { BBMAP_ALIGN                        } from '../../../../modules/nf-core/bbmap/align/main'
include { SELECT_REFERENCE                   } from '../../../../modules/local/select_reference'
include { ASSEMBLY_STRATEGY                  } from '../../../../subworkflows/local/illumina/assembly'


/*
========================================================================================
    SUBWORKFLOW TO PREPROCESSING PIPELINE
========================================================================================
*/

workflow PREPROCESSING_PIPELINE {

    take:
    reads

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    //
    // MODULE: Quality control of raw reads with FastQC
    //
    FASTQC_RAW ( reads )

    //
    // MODULE: Taxonomic classification of raw reads using Kraken2 with PlusPF database
    //
    KRAKEN2_PLUSPF (
        reads,
        params.kraken2_pluspf_database,
        false,
        true
    )

    //
    // MODULE: Adapter trimming with fastp
    //
    FASTP (
        reads,
        params.primer_fasta_RSV_A_and_B,
        false,
        false,
        false
    )

    //
    // MODULE: Quality control of trimmed reads with FastQC
    //
    FASTQC_TRIM ( FASTP.out.reads )

    //
    // MODULE: Remove host reads using Kraken2 with humanDB database
    //
    KRAKEN2_HUMANDB (
        FASTP.out.reads,
        params.kraken2_human_database,
        true,
        true
    )

    //
    // MODULE: Multi-reference alignment with bbmap to select reference
    //
    BBMAP_ALIGN (
        KRAKEN2_HUMANDB.out.unclassified_reads_fastq,
        params.reference_RSV_A_and_B,
        params.bbmap_covstats
    )

    //
    // MODULE: Select best mapping reference
    //
    SELECT_REFERENCE (
        FASTP.out.log,
        BBMAP_ALIGN.out.covstats
    )

    // [meta, fastq] --> [ [id, single_end], fastq1, fastq2, selected_ref ]
    ch_fastq_ref = KRAKEN2_HUMANDB.out.unclassified_reads_fastq
        .join ( SELECT_REFERENCE.out.selected_ref )

    // Channel for RSV-A samples
    // [meta, fastq, RSV-A] --> [ [id, single_end, RSV-A], fastq1, fastq2 ]
    ch_fastq_ref
        .map {
            meta, fastq, selected_ref ->
            if ( selected_ref == 'RSV-A' ) {
                def fmeta = [:]
                fmeta.id = meta.id
                fmeta.single_end = meta.single_end
                fmeta.selected_ref = selected_ref
                [ fmeta, fastq ]
            }
        }
        .set { ch_rsva_pipeline }

    // Channel for RSV-B samples
    // [meta, fastq, RSV-B] --> [ [id, single_end, RSV-B], fastq1, fastq2 ]
    ch_fastq_ref
        .map {
            meta, fastq, selected_ref ->
            if ( selected_ref == 'RSV-B' ) {
                def fmeta = [:]
                fmeta.id = meta.id
                fmeta.single_end = meta.single_end
                fmeta.selected_ref = selected_ref
                [ fmeta, fastq ]
            }
        }
        .set { ch_rsvb_pipeline }

    // Channel for failed MSA samples
    // [meta, fastq, ERROR] --> [ [id, single_end, ERROR], fastq1, fastq2 ]
    ch_fastq_ref
        .map {
            meta, fastq, selected_ref ->
            if ( selected_ref == 'ERROR' ) {
                def fmeta = [:]
                fmeta.id = meta.id
                fmeta.single_end = meta.single_end
                fmeta.selected_ref = selected_ref
                [ fmeta, fastq ]
            }
        }
        .set { ch_assembly_pipeline }

    //
    // SUBWORKFLOW: Assembly strategy subworkflow if previous MSA step failed
    //
    ASSEMBLY_STRATEGY ( ch_assembly_pipeline )

    // Channel for RSV-A samples from assembly strategy
    // [meta, fastq, RSV-A] --> [ [id, single_end, RSV-A], fastq1, fastq2 ]
    ASSEMBLY_STRATEGY.out.fastq_ref
        .map {
            meta, fastq, selected_ref ->
            if ( selected_ref == 'RSV-A' ) {
                def fmeta = [:]
                fmeta.id = meta.id
                fmeta.single_end = meta.single_end
                fmeta.selected_ref = selected_ref
                [ fmeta, fastq ]
            }
        }
        .set { ch_assembly_pipeline_rsva }

    // Channel for RSV-B samples from assembly strategy
    // [meta, fastq, RSV-B] --> [ [id, single_end, RSV-B], fastq1, fastq2 ]
    ASSEMBLY_STRATEGY.out.fastq_ref
        .map {
            meta, fastq, selected_ref ->
            if ( selected_ref == 'RSV-B' ) {
                def fmeta = [:]
                fmeta.id = meta.id
                fmeta.single_end = meta.single_end
                fmeta.selected_ref = selected_ref
                [ fmeta, fastq ]
            }
        }
        .set { ch_assembly_pipeline_rsvb }

    ch_rsva_pipeline
        .concat ( ch_assembly_pipeline_rsva )
        .set { ch_rsva_pipeline_all }

    ch_rsvb_pipeline
        .concat ( ch_assembly_pipeline_rsvb )
        .set { ch_rsvb_pipeline_all }

    // Collect files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_RAW.out.zip.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_PLUSPF.out.report.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP.out.json.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP.out.html.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( FASTP.out.log.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_TRIM.out.zip.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_HUMANDB.out.report.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( ASSEMBLY_STRATEGY.out.multiqc_files )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( FASTQC_RAW.out.versions.first() )
    ch_versions = ch_versions.mix( KRAKEN2_PLUSPF.out.versions.first() )
    ch_versions = ch_versions.mix( FASTP.out.versions.first() )
    ch_versions = ch_versions.mix( FASTQC_TRIM.out.versions.first() )
    ch_versions = ch_versions.mix( KRAKEN2_HUMANDB.out.versions.first() )
    ch_versions = ch_versions.mix( BBMAP_ALIGN.out.versions.first() )
    ch_versions = ch_versions.mix( SELECT_REFERENCE.out.versions.first() )
    ch_versions = ch_versions.mix( ASSEMBLY_STRATEGY.out.versions )

    emit:
    rsva_samples  = ch_rsva_pipeline_all
    rsvb_samples  = ch_rsvb_pipeline_all
    multiqc_files = ch_multiqc_files
    versions      = ch_versions

}
