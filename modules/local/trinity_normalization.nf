process TrinityNormalizeReads {
    tag "$meta.id"
    label 'process_Trinity_Normalization'

    conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.15.1--pl5321h146fbdb_3' }"

    input:
    tuple val(meta), path(reads)
     


    output:
    tuple val(meta), path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files
    path "versions.yml"  , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

     // Define the memory requirements. Trinity needs this as an option.
    def avail_mem = 50
    if (!task.memory) {
        log.info '[Trinity] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

//avail_mem=\$(echo '$avail_mem_groovy')
// prefix=\$(echo '$prefix_groovy')
    '''
    # Ensure the samplesheet.tsv file is empty or create it
    > samplesheet.tsv

    # Loop through all fastq.gz and fq.gz files
    for file in *.{fastq,fq}.gz; do
        # Check if the file actually exists
        if [[ ! -e $file ]]; then
            continue
        fi

        echo "Processing file: $file"

            # Skip if file is an R2 or _2 file
        if [[ $file =~ "_R2.fastq.gz"  || $file =~ "_R2.fq.gz" || $file =~ "_2_val_2.fastq.gz" || $file =~ "_2_val_2.fq.gz" || $file =~ "_2.fastq.gz" || $file =~ "_2.fq.gz" ]]; then
            continue
        fi


        paired_file=""
        # Determine if the file is R1, _1_val_1, or single-end and set paired_file accordingly
        if [[ $file =~ "_R1.fastq.gz" || $file =~ "_R1.fq.gz" ]]; then
             id=$(echo $file | sed 's/_R1.*//')
            paired_file="${id}_R2.${file#*.}"
        elif [[ $file =~ "_1_val_1.fastq.gz" || $file =~ "_1_val_1.fq.gz" ]]; then
            id=$(echo $file | sed 's/_1_val_1.*//')
            paired_file="${id}_2_val_2.${file#*.}"
        elif [[ $file =~ "_1.fastq.gz" || $file =~ "_1.fq.gz" ]]; then
            id=$(echo $file | sed 's/_1.*//')
            paired_file="${id}_2.${file#*.}"      
        else
            id=$(basename "$file" | rev | cut -d "." -f 3- | rev)
        fi

        # Get absolute path for the file
        abs_file=$(realpath "$file")

        # If it's a paired-end read, check if paired_file exists and get its absolute path
        if [[ -n $paired_file && -e $paired_file ]]; then
            abs_paired=$(realpath "$paired_file")
            echo -e "$id\t$id\t$abs_file\t$abs_paired" >> samplesheet.tsv
        elif [[ -z $paired_file ]]; then
            # For single-end reads, only output the single file path
            echo -e "$id\t$id\t$abs_file" >> samplesheet.tsv
        fi
    done
      

        

    Trinity \\
        --seqType fq \\
        --samples_file samplesheet.tsv \\
        --max_memory $avail_mem \\
        --output ${prefix}_trinity \\
        --CPU $task.cpus \\
        --normalize_by_read_set \\
        --just_normalize_reads

     #Use fuzzy matching to find all files matching the *.norm.*.fq pattern.
     #For each matched file, compress its contents using gzip -c and save the compressed contents to a new file with the original filename with the .gz suffix appended via Redirect >.
    #The original file remains unchanged.
    for file in *.norm.*.fq; do
    gzip -cf "$file" > "${file}.gz"
done


    # Adding version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: $(echo $(Trinity --version | head -n 1 2>&1) | sed 's/^Trinity version: Trinity-v//' ))
    END_VERSIONS
'''
}
//In terms of "id=$(echo $file ... part", We will get anything before the following String pattern
//_1_val_1.fastq.gz
//_2_val_2.fastq.gz
//_R1.fastq.gz
//_R2.fastq.gz
//_1.fastq.gz
// _2.fastq.gz
//  
// ''  
//     # Ensure the samplesheet.tsv file is empty or create it
//     > samplesheet.tsv

//     # Loop through all fastq.gz and fq.gz files
//     for file in *.{fastq,fq}.gz; do
//         # Check if the file actually exists
//         if [[ ! -e $file ]]; then
//             continue
//         fi

//         echo "Processing file: $file"

//             # Skip if file is an R2 or _2 file
//         if [[ $file =~ "_R2.fastq.gz"  || $file =~ "_R2.fq.gz" || $file =~ "_2_val_2.fastq.gz" || $file =~ "_2_val_2.fq.gz" || $file =~ "_2.fastq.gz" || $file =~ "_2.fq.gz" ]]; then
//             continue
//         fi


//         paired_file=""
//         # Determine if the file is R1, _1_val_1, or single-end and set paired_file accordingly
//         if [[ $file =~ "_R1.fastq.gz" || $file =~ "_R1.fq.gz" ]]; then
//              id=$(echo $file | sed 's/_R1.*//')
//             paired_file="${id}_R2.${file#*.}"
//         elif [[ $file =~ "_1_val_1.fastq.gz" || $file =~ "_1_val_1.fq.gz" ]]; then
//             id=$(echo $file | sed 's/_1_val_1.*//')
//             paired_file="${id}_2_val_2.${file#*.}"
//         elif [[ $file =~ "_1.fastq.gz" || $file =~ "_1.fq.gz" ]]; then
//             id=$(echo $file | sed 's/_1.*//')
//             paired_file="${id}_2.${file#*.}"      
//         else
//             id=$(basename "$file" | rev | cut -d "." -f 3- | rev)
//         fi

//         # Get absolute path for the file
//         abs_file=$(realpath "$file")

//         # If it's a paired-end read, check if paired_file exists and get its absolute path
//         if [[ -n $paired_file && -e $paired_file ]]; then
//             abs_paired=$(realpath "$paired_file")
//             echo -e "$id\t$id\t$abs_file\t$abs_paired" >> samplesheet.tsv
//         elif [[ -z $paired_file ]]; then
//             # For single-end reads, only output the single file path
//             echo -e "$id\t$id\t$abs_file" >> samplesheet.tsv
//         fi
//     done
