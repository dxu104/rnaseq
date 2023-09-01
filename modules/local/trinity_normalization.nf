//params.samples_file = "/Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/samples.txt" // Path to samples.txt


//nextflow run /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq -profile test,docker --outdir /Users/dxu/MDI  -c /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/conf/base.config -resume

process TrinityNormalizeReads {
    tag "$samples_file"

   //label 'process_high'

    // debug true



        conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.13.2--h00214ad_1' }"
input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    script:
    def samples_list = reads.collect { r ->
        def r1 = r[0]
        def r2 = r.size() > 1 ? r[1] : ""
        "${meta.id}\t${meta.id}_${meta.strandedness}\t${r1}\t${r2}".trim()
    }
    def samples_file_content = samples_list.join("\n")
    def samples_file = "samples.txt"
    samples_file_content.toFile(samples_file)

    // input:
    //  path samples_file

    // Outputs
//    output:
//    path("results_trinity/insilico_read_normalization_altogether/single.norm.*.fq"), emit: single_normalized_files
//    path("results_trinity/insilico_read_normalization_altogether/{left,right}.norm.*.fq"), emit: double_normalized_files
output:
path("results_trinity/insilico_read_normalization_altogether/*.norm.*.fq"), emit: normalized_files

    script:
    """
    Trinity \\
        --seqType fq \\
        --samples_file $samples_file \\
        --max_memory 512G --CPU 64 \\
        --output results_trinity \\
        --normalize_by_read_set \\
        --just_normalize_reads
    """
}
// workflow {
//     TrinityNormalizeReads(params.samples_file)

//     // Access the combined output
//     TrinityNormalizeReads.out.normalized_files.flatten().view {
//         "File: ${it.name}"
//     }
//  --max_memory 10G --CPU 4 \\

// }
