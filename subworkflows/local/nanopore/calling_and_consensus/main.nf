//
// Subworkflow with calling and consensus functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MEDAKA_VARIANTS        } from '../../../../modules/local/medaka_variants'
include { MEDAKA_CONSENSUS       } from '../../../../modules/local/medaka_consensus'

/*
========================================================================================
    SUBWORKFLOW TO PREPROCESSING PIPELINE
========================================================================================
*/

workflow CALLING_AND_CONSENSUS {

    take:
    bam_bai
    reference
    reference_index
    primer_bed
    primer_bed_collapsed
    primer_fasta
    genes_GF_bed

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    //ch_model_variant = DownloadMedakaModel (params.medaka_model_consensus_path, params.medaka_model_consensus_name)

    //
    // MODULE: Calling with Medaka
    //
    MEDAKA_VARIANTS (
        bam_bai,
        reference,
        params.medaka_variant_model
    )

    //
    // MODULE: Consensus with Medaka
    //
    MEDAKA_CONSENSUS (
        bam_bai,
        reference,
        params.medaka_consensus_model
    )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( MEDAKA_VARIANTS.out.versions.first() )
    ch_versions = ch_versions.mix( MEDAKA_CONSENSUS.out.versions.first() )

    emit:
    vcf           = MEDAKA_VARIANTS.out.vcf
    fasta         = MEDAKA_CONSENSUS.out.fasta
    multiqc_files = ch_multiqc_files
    versions      = ch_versions

}
