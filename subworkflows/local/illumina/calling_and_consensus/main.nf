//
// Subworkflow with calling and consensus functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IVAR_VARIANTS        } from '../../../../modules/nf-core/ivar/variants/main'
include { TRANSFORM_TSV_TO_VCF } from '../../../../modules/local/transform_tsv_to_vcf'
include { IVAR_CONSENSUS       } from '../../../../modules/nf-core/ivar/consensus/main'

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

    ch_bam = bam_bai
        .map {
            meta, bam, bai ->
            [ meta, bam ]
        }

    //
    // MODULE: Calling with iVar
    //
    IVAR_VARIANTS (
        ch_bam,
        reference,
        reference_index,
        [],
        true
    )

    //
    // MODULE: Transform iVar TSV to VCF
    //
    TRANSFORM_TSV_TO_VCF (
        IVAR_VARIANTS.out.tsv
    )

    //
    // MODULE: Consensus with iVar
    //
    IVAR_CONSENSUS (
        ch_bam,
        reference,
        true
    )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( IVAR_VARIANTS.out.versions.first() )
    ch_versions = ch_versions.mix( TRANSFORM_TSV_TO_VCF.out.versions.first() )
    ch_versions = ch_versions.mix( IVAR_CONSENSUS.out.versions.first() )

    emit:
    vcf           = TRANSFORM_TSV_TO_VCF.out.vcf
    fasta         = IVAR_CONSENSUS.out.fasta
    multiqc_files = ch_multiqc_files
    versions      = ch_versions

}
