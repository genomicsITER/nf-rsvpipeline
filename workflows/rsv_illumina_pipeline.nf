/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESSING_PIPELINE                              } from '../subworkflows/local/illumina/preprocessing'

include { ALIGN as ALIGN_RSVA                                 } from '../subworkflows/local/illumina/align'
include { ALIGN as ALIGN_RSVB                                 } from '../subworkflows/local/illumina/align'

include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_RSVA } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_RSVB } from '../modules/local/plot_mosdepth_regions'

include { CALLING_AND_CONSENSUS as CALLING_AND_CONSENSUS_RSVA } from '../subworkflows/local/illumina/calling_and_consensus'
include { CALLING_AND_CONSENSUS as CALLING_AND_CONSENSUS_RSVB } from '../subworkflows/local/illumina/calling_and_consensus'

include { FASTA_AGGREGATION as FASTA_AGGREGATION_RSVA         } from '../modules/local/fasta_aggregation'
include { FASTA_AGGREGATION as FASTA_AGGREGATION_RSVB         } from '../modules/local/fasta_aggregation'

include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'

include { paramsSummaryMap                                    } from 'plugin/nf-schema'

include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_rsvpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RSV_ILLUMINA_PIPELINE {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    fastq_exp = params.indir + "/*_R{1,2}_001.fastq.gz"

    ch_fastq = Channel
        .fromFilePairs(fastq_exp)
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    //
    // SUBWORKFLOW: Preprocessing subworkflow
    //
    PREPROCESSING_PIPELINE ( ch_fastq )
    ch_multiqc_files = ch_multiqc_files.mix( PREPROCESSING_PIPELINE.out.multiqc_files )
    ch_versions = ch_versions.mix( PREPROCESSING_PIPELINE.out.versions )

    // * ******************************************* *
    // * RSV-A SAMPLES: ALIGN, CALLING AND CONSENSUS *
    // * ******************************************* *

    //
    // SUBWORKFLOW: Align subworkflow (only for RSV-A samples)
    //
    ALIGN_RSVA (
        PREPROCESSING_PIPELINE.out.rsva_samples,
        params.reference_RSV_A,
        params.reference_RSV_A_index,
        params.primer_bed_RSV_A,
        params.primer_bed_RSV_A_collapsed,
        params.primer_fasta_RSV_A,
        params.genes_GF_bed_RSV_A
    )
    ch_multiqc_files = ch_multiqc_files.mix( ALIGN_RSVA.out.multiqc_files )
    ch_versions = ch_versions.mix( ALIGN_RSVA.out.versions )

    //
    // SUBWORKFLOW: Calling and consensus subworkflow (only for RSV-A samples)
    //
    CALLING_AND_CONSENSUS_RSVA (
        ALIGN_RSVA.out.trimmed_bam_bai,
        params.reference_RSV_A,
        params.reference_RSV_A_index,
        params.primer_bed_RSV_A,
        params.primer_bed_RSV_A_collapsed,
        params.primer_fasta_RSV_A,
        params.genes_GF_bed_RSV_A
    )
    ch_multiqc_files = ch_multiqc_files.mix( CALLING_AND_CONSENSUS_RSVA.out.multiqc_files )
    ch_versions = ch_versions.mix( CALLING_AND_CONSENSUS_RSVA.out.versions )

    //
    // MODULE: Plot mosdepth RSV-A amplicon regions
    //
    PLOT_MOSDEPTH_REGIONS_RSVA (
        ALIGN_RSVA.out.cov_regions_bed.collect { it[1] }
    )
    ch_versions = ch_versions.mix( PLOT_MOSDEPTH_REGIONS_RSVA.out.versions )

    //
    // MODULE: Aggregate RSV-A FASTA files
    //
    FASTA_AGGREGATION_RSVA (
        CALLING_AND_CONSENSUS_RSVA.out.fasta.collect { it[1] }
    )
    ch_versions = ch_versions.mix( FASTA_AGGREGATION_RSVA.out.versions )

    // * ******************************************* *
    // * RSV-B SAMPLES: ALIGN, CALLING AND CONSENSUS *
    // * ******************************************* *

    //
    // SUBWORKFLOW: Align subworkflow (only for RSV-B samples)
    //
    ALIGN_RSVB (
        PREPROCESSING_PIPELINE.out.rsvb_samples,
        params.reference_RSV_B,
        params.reference_RSV_B_index,
        params.primer_bed_RSV_B,
        params.primer_bed_RSV_B_collapsed,
        params.primer_fasta_RSV_B,
        params.genes_GF_bed_RSV_B
    )
    ch_multiqc_files = ch_multiqc_files.mix( ALIGN_RSVB.out.multiqc_files )
    ch_versions = ch_versions.mix( ALIGN_RSVB.out.versions )

    //
    // SUBWORKFLOW: Calling and consensus subworkflow (only for RSV-B samples)
    //
    CALLING_AND_CONSENSUS_RSVB (
        ALIGN_RSVB.out.trimmed_bam_bai,
        params.reference_RSV_B,
        params.reference_RSV_A_index,
        params.primer_bed_RSV_B,
        params.primer_bed_RSV_B_collapsed,
        params.primer_fasta_RSV_B,
        params.genes_GF_bed_RSV_B
    )
    ch_multiqc_files = ch_multiqc_files.mix( CALLING_AND_CONSENSUS_RSVB.out.multiqc_files )
    ch_versions = ch_versions.mix( CALLING_AND_CONSENSUS_RSVB.out.versions )

    //
    // MODULE: Plot mosdepth RSV-B amplicon regions
    //
    PLOT_MOSDEPTH_REGIONS_RSVB (
        ALIGN_RSVB.out.cov_regions_bed.collect { it[1] }
    )
    ch_versions = ch_versions.mix( PLOT_MOSDEPTH_REGIONS_RSVB.out.versions )

    //
    // MODULE: Aggregate RSV-B FASTA files
    //
    FASTA_AGGREGATION_RSVB (
        CALLING_AND_CONSENSUS_RSVB.out.fasta.collect { it[1] }
    )
    ch_versions = ch_versions.mix( FASTA_AGGREGATION_RSVB.out.versions )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rsvpipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()

    summary_params                        = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json" )
    ch_workflow_summary                   = Channel.value( paramsSummaryMultiqc(summary_params) )
    ch_multiqc_files                      = ch_multiqc_files.mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value( methodsDescriptionText(ch_multiqc_custom_methods_description) )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix( ch_methods_description.collectFile( name: 'methods_description_mqc.yaml', sort: true ) )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
