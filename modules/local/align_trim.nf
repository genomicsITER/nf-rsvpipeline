process ALIGN_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::artic=1.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.1--pyhdfd78af_0' :
        'quay.io/biocontainers/artic:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def report_csv = "${prefix}.alignreport.csv"
    def report_err = "${prefix}.alignreport.err"
    def trimmed_bam = "${prefix}.sorted.bam"
    """
    align_trim $bed --trim-primers --report $report_csv --verbose < $bam 2> $report_err | \
    samtools sort -T ${meta.id} - -o $trimmed_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$( artic -v | cut -d' ' -f2 )
    END_VERSIONS
    """
}
