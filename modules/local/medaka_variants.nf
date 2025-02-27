process MEDAKA_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::medaka=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.0.1--py38h8774169_0' :
        'quay.io/biocontainers/medaka:2.0.1--py38h8774169_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    val(model)

    output:
    tuple val(meta), path('*.vcf')   , emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    medaka_variant -m $model -i $bam -r $reference

    cp medaka/medaka.annotated.vcf ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | cut -d' ' -f2 )
    END_VERSIONS
    """
}
