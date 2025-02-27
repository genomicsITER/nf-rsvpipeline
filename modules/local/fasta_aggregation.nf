process FASTA_AGGREGATION {
    label 'process_medium'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'quay.io/biocontainers/ubuntu:22.04' }"

    input:
    path fastas

    output:
    path '*_samples.consensus.fa', emit: aggregated_fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "aggregated_fasta"
    """
    cat $fastas > $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version | head -1 | sed 's/cat (GNU coreutils) //g')
    END_VERSIONS
    """
}
