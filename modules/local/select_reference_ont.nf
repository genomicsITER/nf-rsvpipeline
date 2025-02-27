process SELECT_REFERENCE_ONT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), env('selected_ref')           , emit: selected_ref
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    strain=\$( grep '>' $fasta | head -1 | sed 's/>//g' )

    if [[ \${strain} == *"RSV_A"* || \${strain} == *"RSV_AD"* ]]; then
        selected_ref=\$( echo "RSV-A" )
    elif [[ \${strain} == *"RSV_B"* || \${strain} == *"RSV_BD"* ]]; then
        selected_ref=\$( echo "RSV-B" )
    else
        selected_ref=\$( echo "ERROR" )
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
