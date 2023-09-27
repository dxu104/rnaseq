process TRINITY_NORMALIZATION {
    tag "$meta.id"
    label 'process_Trinity_Normalization_Parallel'

    conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.15.1--pl5321h146fbdb_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trinity/*/*.fq.gz")       , emit: transcript_fastq
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        reads_args = "--single ${reads}"
    } else {
        reads_args = "--left ${reads[0]} --right ${reads[1]}"
    }

    // --seqType argument, fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    seqType_args = reads[0] ==~ /(.*fasta(.gz)?$)|(.*fa(.gz)?$)/ ? "fa" : "fq"

    // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = 50
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.95).intValue()
    }

    """
    # Note that Trinity needs the word 'trinity' in the outdir

    Trinity \\
    --seqType ${seqType_args} \\
    --max_memory ${avail_mem}G \\
    ${reads_args} \\
    --output ${prefix}_trinity \\
    --CPU $task.cpus \\
    --normalize_max_read_cov 30           \\
    --just_normalize_reads\\
    $args

    #this is from offical template, do not use gzip -cf ${prefix}_trinity.Trinity.fasta > ${prefix}.fa.gz

    #Use fuzzy matching to find all files matching the *.fq pattern.
    #For each matched file, compress its contents using gzip -c and save the compressed contents to a new file with the original filename with the .gz suffix appended via Redirect >.
    #The original file remains unchanged.
    # Check for files with the desired pattern and gzip them
    # Check for files with the desired pattern in current directory and subdirectories, then gzip them
    found_files=\$(find . -type f -name "*.fq")
    if [[ -n "\$found_files" ]]; then
        find . -type f -name "*.fq" | while read -r file; do
            echo "Processing file: \$file"
            gzip -cf "\$file" > "\${file}.gz"
            echo "Compressed to: \${file}.gz"
        done
    else
        echo "No files matching the pattern were found."
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(echo \$(Trinity --version | head -n 1 2>&1) | sed 's/^Trinity version: Trinity-v//' ))
    END_VERSIONS

    # Need to only take the first line of --version since it will warn about not being up-to-date and this messes up the version.yaml.
    """
}
