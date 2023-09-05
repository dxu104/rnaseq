//params.samples_file = "/Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/samples.txt" // Path to samples.txt


//nextflow run /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq -profile test,docker --outdir /Users/dxu/MDI  -c /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/conf/base.config -resume


process TrinityNormalizeReads_SingleEnd {

    conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.15.1--pl5321h146fbdb_3' }"

    input:
    tuple path(read1)

    output:
    path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files



    script:
    '''
# Create an array of file paths for R1 files
r1_files=(*_R1.fastq.gz)

# Ensure the samplesheet.tsv file is empty or create it
> samplesheet.tsv

# Loop through the R1 files and find the matching R2 file
for r1 in ${r1_files[@]}; do
    # Extract the ID from the filename
    id=$(echo $r1 | cut -f 1 -d "_")

    # Build the corresponding R2 filename
    r2="${id}_R2.fastq.gz"

    # Get absolute paths for each file
    abs_r1=$(realpath "$r1")
    abs_r2=$(realpath "$r2")

    # Append information to the TSV file
    echo -e "$id\t$id\t$abs_r1\t$abs_r2" >> samplesheet.tsv
done

Trinity \\
    --seqType fq \\
    --samples_file samplesheet.tsv \\
    --max_memory 512G --CPU 64 \\
    --output results_trinity \\
    --normalize_by_read_set \\
    --just_normalize_reads
'''
}






// def read1_files, read2_files
    // if (meta.single_end == true) {
    //     read1_files = reads.collect { it.toString() }.join(",")
    //     read2_files = null
    // } else {
    //     read1_files = reads[0].collect { it.toString() }.join(",")
    //     read2_files = reads[1].collect { it.toString() }.join(",")
    // }

// """
    // if [ "${meta.single_end}" == "true" ]
    // then
    //     Trinity \\
    //         --seqType fq \\
    //         --single $read1_files \\
    //         --max_memory 512G --CPU 32 \\
    //         --output results_trinity \\
    //         --normalize_by_read_set \\
    //         --just_normalize_reads
    // else
    //     echo "Read2 Files: $read2_files"
    //     Trinity \\
    //         --seqType fq \\
    //         --left $read1_files \\
    //         --right $read2_files \\
    //         --max_memory 512G --CPU 32 \\
    //         --output results_trinity \\
    //         --normalize_by_read_set \\
    //         --just_normalize_reads
    // fi
    // """

// 下面这个生效了，但是生成了7个left,7个right。
// process TrinityNormalizeReads {
//     //tag "$samples_file"

//    //label 'process_high'

//    // debug true

//     conda "bioconda::trinity=2.13.2"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
//         'biocontainers/trinity:2.13.2--h00214ad_1' }"

//           input:
//     tuple val(meta), path(reads)

//         output:
//     path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files

//    script:
// def read1 = reads[0].toString()
// def read2 = reads[1].toString()
// """
// echo "Sample ID: ${meta.id}"
// echo "Read1 Path: $read1"
// echo "Read2 Path: $read2"
//     Trinity \\
//         --seqType fq \\
//          --left $read1 \\
//          --right $read2 \\
//         --max_memory 512G --CPU 64 \\
//         --output results_trinity \\
//         --normalize_by_read_set \\
//         --just_normalize_reads

// """


// }

// input:
//     tuple val(meta), path(reads)
// //, stageAs: "input*/*"
//     output:
//     //path("samples_double_end.txt"), emit: samples
//     //path("samples.txt"), emit: file_samples
// path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files

//     script:
// ch_samples_file=Channel.fromPath(reads).map { path ->
//         def read1 = path[0].toString()
//         def read2 = path.size() > 1 ? path[1].toString() : ""
//         return [meta.id, "${meta.id}_${meta.strandedness}", read1, read2].join('\t')
//     }
//     .collectFile(name: 'samples_file.txt', newLine: true, sort: true)


//     ch_samples_file.view { txt ->
//     "Sample File Content: $txt"
// }

//     """
    // Trinity \\
    //     --seqType fq \\
    //     --samples_file $ch_samples_file\\
    //     --max_memory 512G --CPU 64 \\
    //     --output results_trinity \\
    //     --normalize_by_read_set \\
    //     --just_normalize_reads
    // """
// }







//     input:
//      path samples_file

//     Outputs
//    output:
//    path("results_trinity/insilico_read_normalization_altogether/single.norm.*.fq"), emit: single_normalized_files
//    path("results_trinity/insilico_read_normalization_altogether/{left,right}.norm.*.fq"), emit: double_normalized_files

// process TrinityNormalizeReads {
//     conda "bioconda::trinity=2.13.2"
//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
//         'biocontainers/trinity:2.13.2--h00214ad_1' }"

//     input:
//     path samples_file
//     output:
//     path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files

//     script:
//     """
//     Trinity \\
//         --seqType fq \\
//         --samples_file $samples_file \\
//         --max_memory 512G --CPU 64 \\
//         --output results_trinity \\
//         --normalize_by_read_set \\
//         --just_normalize_reads
//     """

