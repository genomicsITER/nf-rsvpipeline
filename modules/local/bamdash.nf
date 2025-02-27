process BAMDASH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bamdash=0.3.1-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamdash:0.3.1--pyhdfd78af_0' :
        'quay.io/biocontainers/bamdash:0.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bam_header)
    path(bed)
    path(bed2)

    output:
    tuple val(meta), path('*.html') , emit: html
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ref_id=\$( grep "@SQ" $bam_header | awk '{ print \$2 }' | cut -d":" -f2 )

    bamdash -b $bam -r "\${ref_id}" -t $bed $bed2

    mv plot.html ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamdash: \$( bamdash --version | cut -d' ' -f2 )
    END_VERSIONS
    """
}
