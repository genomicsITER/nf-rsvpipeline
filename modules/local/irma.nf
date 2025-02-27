process IRMA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::irma=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/irma:1.2.0--pl5321hdfd78af_0' :
        'quay.io/biocontainers/irma:1.2.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    val(module)

    output:
    tuple val(meta), path('*.fasta') , emit: fasta
    tuple val(meta), path('*.bam')   , emit: bam
    tuple val(meta), path('*.bai')   , emit: bai
    tuple val(meta), path('*.vcf')   , emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    IRMA \
        $module \
        $fastq \
        $prefix/${meta.id}

    mv $prefix/${meta.id}/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irma: \$( IRMA | head -1 | cut -d',' -f2 | cut -d' ' -f2 | sed 's/v//g' )
    END_VERSIONS
    """
}
