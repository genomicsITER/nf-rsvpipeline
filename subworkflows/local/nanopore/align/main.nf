//
// Subworkflow with align functionality specific to the rsvpipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN                                     } from '../../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM_ONT           } from '../../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_ONT         } from '../../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_ONT         } from '../../../../modules/nf-core/samtools/idxstats/main'

include { SAMTOOLS_VIEW as DISCARD_UNMAPPED_READS_ONT        } from '../../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MAPPED_BAM_ONT    } from '../../../../modules/nf-core/samtools/index/main'

include { ALIGN_TRIM                                         } from '../../../../modules/local/align_trim'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRIMMED_BAM_ONT   } from '../../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TRIMMED_ONT } from '../../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_TRIMMED_ONT } from '../../../../modules/nf-core/samtools/idxstats/main'

include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_ONT               } from '../../../../modules/nf-core/samtools/depth/main'
include { MOSDEPTH as MOSDEPTH_ONT                           } from '../../../../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_PER_AMPLICON_ONT              } from '../../../../modules/nf-core/mosdepth/main'

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
    // MODULE: Align reads with Minimap2
    //
    MINIMAP2_ALIGN (
        reads,
        ch_reference_files.map {
            meta, reference, primer_bed, primer_bed_collapsed, primer_fasta, genes_GF_bed ->
            [meta, reference]
        },
        true,
        "bai",
        false,
        false
    )

    // [meta, bam, bai]
    ch_bam_bai = MINIMAP2_ALIGN.out.bam
        .join ( MINIMAP2_ALIGN.out.index )

    //
    // MODULES: Quality control with Samtools flagstat and idxstats
    //
    SAMTOOLS_FLAGSTAT_ONT ( ch_bam_bai )
    SAMTOOLS_IDXSTATS_ONT ( ch_bam_bai )

    //
    // MODULE: Discard unmapped reads with Samtools view
    //
    DISCARD_UNMAPPED_READS_ONT (
        ch_bam_bai,
        [[], []],
        []
    )

    //
    // MODULE: Index sorted and mapped bam
    //
    SAMTOOLS_INDEX_MAPPED_BAM_ONT ( DISCARD_UNMAPPED_READS_ONT.out.bam )

    // [meta, mapped_bam, bai]
    ch_mapped_bam_bai = DISCARD_UNMAPPED_READS_ONT.out.bam
        .join ( SAMTOOLS_INDEX_MAPPED_BAM_ONT.out.bai )

    //
    // MODULE: Trim adapters with align_trim (from ARTIC) before create consensus sequences
    //
    ALIGN_TRIM (
        ch_mapped_bam_bai,
        primer_bed
    )

    //
    // MODULE: Index trimmed bam
    //
    SAMTOOLS_INDEX_TRIMMED_BAM_ONT ( ALIGN_TRIM.out.bam )

    // [meta, bam, bai]
    ch_trimmed_bam_bai = ALIGN_TRIM.out.bam
        .join ( SAMTOOLS_INDEX_TRIMMED_BAM_ONT.out.bai )

    //
    // MODULES: Quality control with Samtools flagstat, idxstats and depth
    //
    SAMTOOLS_FLAGSTAT_TRIMMED_ONT ( ch_trimmed_bam_bai )
    SAMTOOLS_IDXSTATS_TRIMMED_ONT ( ch_trimmed_bam_bai )
    SAMTOOLS_DEPTH_ONT ( ALIGN_TRIM.out.bam, [[], []] )

    // [meta, bam, bai, null]
    ch_mosdepth = ch_trimmed_bam_bai
        .map {
            meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }

    //
    // MODULE: Coverage analysis with Mosdepth (all positions)
    //
    MOSDEPTH_ONT (
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
    MOSDEPTH_PER_AMPLICON_ONT (
        ch_mosdepth_amplicon,
        [[], []]
    )

    // Collect files for MultiQC
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_ONT.out.flagstat.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_IDXSTATS_ONT.out.idxstats.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_TRIMMED_ONT.out.flagstat.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_IDXSTATS_TRIMMED_ONT.out.idxstats.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_DEPTH_ONT.out.tsv.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_ONT.out.global_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_ONT.out.summary_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_ONT.out.regions_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON_ONT.out.global_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON_ONT.out.summary_txt.collect{it[1]} )
    ch_multiqc_files = ch_multiqc_files.mix( MOSDEPTH_PER_AMPLICON_ONT.out.regions_txt.collect{it[1]} )

    // Collect versions for MultiQC
    ch_versions = ch_versions.mix( MINIMAP2_ALIGN.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_IDXSTATS_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( DISCARD_UNMAPPED_READS_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_MAPPED_BAM_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( ALIGN_TRIM.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_TRIMMED_BAM_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT_TRIMMED_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_IDXSTATS_TRIMMED_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( SAMTOOLS_DEPTH_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( MOSDEPTH_ONT.out.versions.first() )
    ch_versions = ch_versions.mix( MOSDEPTH_PER_AMPLICON_ONT.out.versions.first() )

    emit:
    trimmed_bam_bai = ch_trimmed_bam_bai
    cov_regions_bed = MOSDEPTH_PER_AMPLICON_ONT.out.regions_bed
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions

}
