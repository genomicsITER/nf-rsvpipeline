//
// Subworkflow with assembly functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPADES                  } from '../../../../modules/nf-core/spades/main'
include { BLAST_BLASTN            } from '../../../../modules/nf-core/blast/blastn/main'
include { SELECT_REFERENCE_BLASTN } from '../../../../modules/local/select_reference_blastn'


/*
========================================================================================
    SUBWORKFLOW TO PREPROCESSING PIPELINE
========================================================================================
*/

workflow ASSEMBLY_STRATEGY {

    take:
    reads

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    ch_spades = reads
        .map {
            meta, reads -> [ meta, reads, [], [] ]
        }

    //
    // MODULE: Assembly with Spades
    //
    SPADES (
        ch_spades,
        [],
        []
    )

    //
    // MODULE: Detect hits using BLASTn
    //
    BLAST_BLASTN (
        SPADES.out.scaffolds,
        params.blast_database
    )

    //
    // MODULE: Select reference based on blastn results
    //
    SELECT_REFERENCE_BLASTN (
        BLAST_BLASTN.out.txt
    )

    // [meta, fastq] --> [ [id, single_end, ERROR], fastq1, fastq2, selected_ref ]
    ch_fastq_ref = reads
        .join ( SELECT_REFERENCE_BLASTN.out.selected_ref )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( SPADES.out.versions.first() )
    ch_versions = ch_versions.mix( BLAST_BLASTN.out.versions.first() )
    ch_versions = ch_versions.mix( SELECT_REFERENCE_BLASTN.out.versions.first() )

    emit:
    fastq_ref     = ch_fastq_ref
    multiqc_files = ch_multiqc_files
    versions      = ch_versions

}
