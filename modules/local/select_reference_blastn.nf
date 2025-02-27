process SELECT_REFERENCE_BLASTN {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(blastn_txt)

    output:
    tuple val(meta), env('selected_ref')           , emit: selected_ref
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Parse BLASTn results:
    selected_reference=\$( sort -t\$'\\t' -k3 -k4 -nr $blastn_txt | head -1 )

    if [[ \${selected_reference} == *"Human respiratory syncytial virus A"* ]]; then
        selected_ref=\$( echo "RSV-A" )
    elif [[ \${selected_reference} == *"Human respiratory syncytial virus B"* ]]; then
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
