//
// Subworkflow with align functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX                                      } from '../../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                                        } from '../../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM           } from '../../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT                              } from '../../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS                              } from '../../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_VIEW as DISCARD_UNMAPPED_READS        } from '../../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MAPPED_BAM    } from '../../../../modules/nf-core/samtools/index/main'

include { IVAR_TRIM                                      } from '../../../../modules/nf-core/ivar/trim/main'
include { SAMTOOLS_SORT                                  } from '../../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRIMMED_BAM   } from '../../../../modules/nf-core/samtools/index/main'

include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TRIMMED } from '../../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_TRIMMED } from '../../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH                                 } from '../../../../modules/nf-core/samtools/depth/main'

include { MOSDEPTH                                       } from '../../../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_PER_AMPLICON              } from '../../../../modules/nf-core/mosdepth/main'

/*
========================================================================================
    SUBWORKFLOW TO PREPROCESSING PIPELINE
========================================================================================
*/

workflow ALIGN {

    take:
    reads
    reference
    reference_index
    primer_bed
    primer_bed_collapsed
    primer_fasta
    genes_GF_bed

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    // [meta, reads] --> [ [id, single_end, selected_ref], reference, primer_bed, primer_bed_collapsed, primer_fasta, genes_GF_bed ]
    ch_reference_files = reads
        .map {
            meta, reads ->
            [ meta, reference, primer_bed, primer_bed_collapsed, primer_fasta, genes_GF_bed ]
        }

    //
    // MODULE: Index reference
    //
    BWA_INDEX (
        ch_reference_files.map {
            meta, reference, primer_bed, primer_bed_collapsed, primer_fasta, genes_GF_bed ->
            [meta, reference]
        }
    )

    //
    // MODULE: Align reads with BWA_MEM
    //
    BWA_MEM (
        reads,
        BWA_INDEX.out.index,
        ch_reference_files.map {
            meta, reference, primer_bed, primer_bed_collapsed, primer_fasta, genes_GF_bed ->
            [meta, reference]
        },
        true
    )

    //
    // MODULE: Index sorted bam
    //
    SAMTOOLS_INDEX_BAM ( BWA_MEM.out.bam )

    // [meta, bam, bai]
    ch_bam_bai = BWA_MEM.out.bam
        .join ( SAMTOOLS_INDEX_BAM.out.bai )

    //
    // MODULES: Quality control with Samtools flagstat and idxstats
    //
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    SAMTOOLS_IDXSTATS ( ch_bam_bai )

    //
    // MODULE: Discard unmapped reads with Samtools view
    //
    DISCARD_UNMAPPED_READS (
        ch_bam_bai,
        [[], []],
        []
    )

    //
    // MODULE: Index sorted and mapped bam
    //
    SAMTOOLS_INDEX_MAPPED_BAM ( DISCARD_UNMAPPED_READS.out.bam )

    // [meta, mapped_bam, bai]
    ch_mapped_bam_bai = DISCARD_UNMAPPED_READS.out.bam
        .join ( SAMTOOLS_INDEX_MAPPED_BAM.out.bai )

    //
    // MODULE: Trim adapters with iVar before create consensus sequences
    //
    IVAR_TRIM (
        ch_mapped_bam_bai,
        params.primer_bed_RSV_A_and_B
    )

    //
    // MODULE: Sort trimmed bam
    //
    SAMTOOLS_SORT (
        IVAR_TRIM.out.bam,
        [[], []]
    )

    //
    // MODULE: Index trimmed bam
    //
    SAMTOOLS_INDEX_TRIMMED_BAM ( SAMTOOLS_SORT.out.bam )

    // [meta, bam, bai]
    ch_trimmed_bam_bai = SAMTOOLS_SORT.out.bam
        .join ( SAMTOOLS_INDEX_TRIMMED_BAM.out.bai )

    //
    // MODULES: Quality control with Samtools flagstat, idxstats and depth
    //
    SAMTOOLS_FLAGSTAT_TRIMMED ( ch_trimmed_bam_bai )
    SAMTOOLS_IDXSTATS_TRIMMED ( ch_trimmed_bam_bai )
    SAMTOOLS_DEPTH ( SAMTOOLS_SORT.out.bam, [[], []] )

    // [meta, bam, bai, null]
    ch_mosdepth = ch_trimmed_bam_bai
        .map {
            meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }

    //
    // MODULE: Coverage analysis with Mosdepth (all positions)
    //
    MOSDEPTH (
        ch_mosdepth,
        [[], []]
    )

    // [meta, bam, bai, bed]
    ch_mosdepth_amplicon = ch_trimmed_bam_bai
        .map {
            meta, bam, bai ->
            [ meta, bam, bai, primer_bed_collapsed ]
        }

    //
    // MODULE: Coverage analysis with Mosdepth (amplicon regions)
    //
    MOSDEPTH_PER_AMPLICON (
        ch_mosdepth_amplicon,
        [[], []]
    )

    // Collect files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_TRIMMED.out.flagstat.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_IDXSTATS_TRIMMED.out.idxstats.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_DEPTH.out.tsv.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH.out.global_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH.out.summary_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH.out.regions_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON.out.global_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON.out.summary_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON.out.regions_txt.collect{it[1]} )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( BWA_INDEX.out.versions.first() )
    ch_versions = ch_versions.mix( BWA_MEM.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_BAM.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_IDXSTATS.out.versions.first() )
    ch_versions = ch_versions.mix( DISCARD_UNMAPPED_READS.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_MAPPED_BAM.out.versions.first() )
    ch_versions = ch_versions.mix( IVAR_TRIM.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_TRIMMED_BAM.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT_TRIMMED.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_IDXSTATS_TRIMMED.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH.out.versions.first() )
    ch_versions = ch_versions.mix( MOSDEPTH.out.versions.first() )
    ch_versions = ch_versions.mix( MOSDEPTH_PER_AMPLICON.out.versions.first() )

    emit:
    trimmed_bam_bai = ch_trimmed_bam_bai
    cov_regions_bed = MOSDEPTH_PER_AMPLICON.out.regions_bed
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions

}
