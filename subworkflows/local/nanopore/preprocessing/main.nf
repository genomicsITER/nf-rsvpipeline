//
// Subworkflow with preprocessing functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NANOPLOT as NANOPLOT_RAW               } from '../../../../modules/nf-core/nanoplot/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_PLUSPF_ONT  } from '../../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HUMANDB_ONT } from '../../../../modules/nf-core/kraken2/kraken2/main'
include { NANOPLOT as NANOPLOT_SCRUBBED          } from '../../../../modules/nf-core/nanoplot/main'
include { IRMA                                   } from '../../../../modules/local/irma'
include { SELECT_REFERENCE_ONT                   } from '../../../../modules/local/select_reference_ont'

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
    // MODULE: Quality control of raw reads with NanoPlot
    //
    NANOPLOT_RAW ( reads )

    //
    // MODULE: Taxonomic classification of raw reads using Kraken2 with PlusPF database
    //
    KRAKEN2_PLUSPF_ONT (
        reads,
        params.kraken2_pluspf_database,
        false,
        true
    )

    //
    // MODULE: Remove host reads using Kraken2 with humanDB database
    //
    KRAKEN2_HUMANDB_ONT (
        reads,
        params.kraken2_human_database,
        true,
        true
    )

    //
    // MODULE: Quality control of scrubbed reads with NanoPlot
    //
    NANOPLOT_SCRUBBED ( KRAKEN2_HUMANDB_ONT.out.unclassified_reads_fastq, )

    //
    // MODULE: Strain assignment using IRMA
    //
    IRMA (
        KRAKEN2_HUMANDB_ONT.out.unclassified_reads_fastq,
        params.irma_module
    )

    //
    // MODULE: Select best mapping reference
    //
    SELECT_REFERENCE_ONT (
        IRMA.out.fasta
    )

    // [meta, fastq] --> [ [id, single_end], fastq, selected_ref ]
    ch_fastq_ref = KRAKEN2_HUMANDB_ONT.out.unclassified_reads_fastq
        .join ( SELECT_REFERENCE_ONT.out.selected_ref )

    // Channel for RSV-A samples
    // [meta, fastq, RSV-A] --> [ [id, single_end, RSV-A], fastq ]
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
    // [meta, fastq, RSV-B] --> [ [id, single_end, RSV-B], fastq ]
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

    // Collect files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix( NANOPLOT_RAW.out.txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_PLUSPF_ONT.out.report.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_HUMANDB_ONT.out.report.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( NANOPLOT_SCRUBBED.out.txt.collect{it[1]} )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( NANOPLOT_RAW.out.versions.first() )
    ch_versions = ch_versions.mix( KRAKEN2_PLUSPF_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( KRAKEN2_HUMANDB_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( NANOPLOT_SCRUBBED.out.versions.first() )
    ch_versions = ch_versions.mix( IRMA.out.versions.first() )
    ch_versions = ch_versions.mix( SELECT_REFERENCE_ONT.out.versions.first() )

    emit:
    rsva_samples  = ch_rsva_pipeline
    rsvb_samples  = ch_rsvb_pipeline
    multiqc_files = ch_multiqc_files
    versions      = ch_versions

}
