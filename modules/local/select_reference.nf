process SELECT_REFERENCE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fastp_log)
    tuple val(meta), path(covstats)

    output:
    tuple val(meta), env('selected_ref')           , emit: selected_ref
    tuple val(meta), path('*_refs.tsv')            , optional:true, emit: refs
    tuple val(meta), path('*_failed_assembly.tsv') , optional:true, emit: failed
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    raw_reads=\$( grep -A1 "before filtering:" $fastp_log | grep "total reads:" | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}' )
    trimmed_reads=\$( grep -A1 "after filtering:" $fastp_log | grep "total reads:" | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}' )

    # Select best mapping:
    m=3
    p=0

    # Select reference using this script from Graniger-Lab's Revica pipeline:
    # https://github.com/greninger-lab/revica/blob/main/bin/select_reference.py

    select_reference.py \
        -bbmap_covstats $covstats \
        -raw_reads \${raw_reads} \
        -trimmed_reads \${trimmed_reads} \
        -b ${meta.id} \
        -m \${m} \
        -p \${p}

    if [ -f ${meta.id}_refs.tsv ]; then
        # Parse previous output file to select final reference to use in the following steps:
        selected_reference=\$( cat ${meta.id}_refs.tsv )

        if [[ \${selected_reference} == *"hRSV_A"* ]]; then
            selected_ref=\$( echo "RSV-A" )
        elif [[ \${selected_reference} == *"hRSV_B"* ]]; then
            selected_ref=\$( echo "RSV-B" )
        else
            selected_ref=\$( echo "ERROR" )
        fi
    else
        selected_ref=\$( echo "ERROR" )
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
