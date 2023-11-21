/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['star_salmon', 'star_rsem', 'hisat2'],
    trimmers       : ['trimgalore', 'fastp'],
    pseudoaligners : ['salmon'],
    rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.ribo_database_manifest, params.splicesites,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if (params.bam_csi_index) {
    for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
        if (rseqc_modules.contains(rseqc_module)) {
            rseqc_modules.remove(rseqc_module)
        }
    }
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc'
include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/local/umitools_prepareforrsem.nf'
//customized modules
// include { TrinityNormalizeReads as TrinityNormalizeReads_SingleEnd } from '../modules/local/TrinityNormalization/trinity_normalization.nf'
// include { TrinityNormalizeReads as TrinityNormalizeReads_DoubleEnd } from '../modules/local/TrinityNormalization/trinity_normalization.nf'
//include { CreateSampleFile } from '../modules/local/Samples_file_for_trinity_normalization.nf'
// include { Staging as Staging_SingleEnd } from '../modules/local/create_samples_file_staging.nf'
// include { Staging as Staging_DoubleEnd } from '../modules/local/create_samples_file_staging.nf'

include { BAMSIFTER  } from '../modules/local/bamsifter.nf'
include { BAMSIFTER as BAMSIFTER_NORMALIZATION_MERGED_BAM} from '../modules/local/bamsifter.nf'
include { SAMTOOLS_MERGE } from '../modules/local/samtools_merge.nf'
// include { TRINITY_NORMALIZATION as TRINITY_NORMALIZATION_PARALLEL_DoubleEnd} from '../modules/local/trinity.nf'
// include { TRINITY_NORMALIZATION as TRINITY_NORMALIZATION_PARALLEL_SingleEnd} from '../modules/local/trinity.nf'

include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA_FROM_NEW_GTF      } from '../modules/nf-core/rsem/preparereference/main'
include { GTF_GENE_FILTER   as   GTF_GENE_FILTER_FROM_NEW_GTF                    } from '../modules/local/gtf_gene_filter.nf'
include { SALMON_INDEX   as   SALMON_INDEX_FROM_NEW_TRANSCRIPT_FASTA                 } from '../modules/nf-core/salmon/index/main'
include { GFFCOMPARE                   } from '../modules/local/gffcompare.nf'

//fastq after trinity normalization
//include { FASTQC as FASTQC_AFTER_TRINITY} from '../modules/nf-core/fastqc/main'

//StringTie merge modules
include {STRINGTIE_MERGE} from '../modules/local/stringTie_merge/main.nf'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { TRINITY_NORMALIZATION  } from '../subworkflows/local/trinity_normalization'
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR     } from '../subworkflows/local/align_star'
include { QUANTIFY_RSEM  } from '../subworkflows/local/quantify_rsem'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'

//!! we add this
include { SAMTOOLS_SORT as   SAMTOOLS_SORT_BAM  } from '../modules/nf-core/samtools/sort/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ             } from '../modules/nf-core/qualimap/rnaseq/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { STRINGTIE_STRINGTIE         } from '../modules/nf-core/stringtie/stringtie/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE as FASTQ_FASTQC_UMITOOLS_TRIMGALORE_AFTER_TRINITY } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_ALIGN_HISAT2               } from '../subworkflows/nf-core/fastq_align_hisat2/main'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_RSEQC                        } from '../subworkflows/nf-core/bam_rseqc/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow RNASEQ_TRANSCRIPTOME_UPDATE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.gencode,
        is_aws_igenome,
        biotype,
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { WorkflowRnaseq.checkMaxContigSize(it, log) }
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    PREPARE_GENOME.out.fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        !params.salmon_index && !('salmon' in prepareToolIndices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: WorkflowRnaseq.getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    // default setting trimmer is params.trimmer = 'trimgalore' in the nextflow.config
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            PREPARE_GENOME.out.bbsplit_index,
            [],
            [ [], [] ],
            false
        )
        .primary_fastq
        .set { ch_filtered_reads }
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
//why we can use collect() to get ch_sortmerna_fastas, and pass ch_sortmerna_fastas into SORTMERNA?
// but do not use collect() for ch_filtered_reads because collect will broke the tuple stucture
//and collect() will not broke the structure for ch_sortmerna_fastas
        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .view{ "Meta: ${it[0]}, Path: ${it[1]}" }
        .set { ch_filtered_reads }

        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
}


/*

Channel.empty().set { ch_normalized_double_end_files }

if (params.double_end_sample) {

ch_samples_double_end = ch_filtered_reads
    .filter { meta, path ->
        meta.single_end == false || meta.single_end == null  // null is for the case of undefined
    }

     //simultaneously run Trinity to expedite the run time
    TRINITY_NORMALIZATION_PARALLEL_DoubleEnd(ch_samples_double_end)

    // collect() to collect all  transcript_fasta
    ch_samples_double_end= TRINITY_NORMALIZATION_PARALLEL_DoubleEnd.out.transcript_fastq.collect()

//

ch_versions = ch_versions.mix(TRINITY_NORMALIZATION_PARALLEL_DoubleEnd.out.versions)

//After collect Double End Sample Text Content: [[id:RAP1_IAA_30M_REP1, single_end:false, strandedness:reverse], [/mdibl-nextflow-work/dxu/smallestTest_09-21-23_memvergeOndemand/53/b37e086f0061c70d750d1b453eac97/RAP1_IAA_30M_REP1_trinity/insilico_read_normalization/RAP1_IAA_30M_REP1_1.non_rRNA.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq.gz, /mdibl-nextflow-work/dxu/smallestTest_09-21-23_memvergeOndemand/53/b37e086f0061c70d750d1b453eac97/RAP1_IAA_30M_REP1_trinity/insilico_read_normalization/RAP1_IAA_30M_REP1_2.non_rRNA.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq.gz], [id:WT_REP2, single_end:false, strandedness:reverse], [/mdibl-nextflow-work/dxu/smallestTest_09-21-23_memvergeOndemand/7a/c567da50264ddcc70e57c27720bdb7/WT_REP2_trinity/insilico_read_normalization/WT_REP2_1.non_rRNA.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq.gz, /mdibl-nextflow-work/dxu/smallestTest_09-21-23_memvergeOndemand/7a/c567da50264ddcc70e57c27720bdb7/WT_REP2_tMonitor the execution with Nextflow Tower using this URL: https://tower.nf/orgs/MDIBL-Biocore/workspaces/Memverge/watch/5R3w7digPTZxbMexecutor >  float (32)[ac/a3170f] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔[e2/d4240c] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔[2a/343603] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔[a9/434460] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔[7b/5cc81d] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔[d0/6cce9c] process > NFCORE_RNASEQ:RN... [100%] 1 of 1 ✔
            ch_samples_double_end
            .flatMap()
            .map {
                 item ->
            if(item instanceof Map) {
            item['id'] = 'all_double'
            return item
            } else {
            return item
            }
            }
            .collate(2)
            .groupTuple()
            .map {
            meta, fastq ->
            [ meta, fastq.flatten() ]
            }
            .set{ch_inputfor_double_TrinityNormalization}
            ch_inputfor_double_TrinityNormalization.view{ "Ready for Second Trinity  Meta: ${it[0]}, Path: ${it[1]}" }




TrinityNormalizeReads_DoubleEnd(ch_inputfor_double_TrinityNormalization)


 ch_normalized_double_end_files = TrinityNormalizeReads_DoubleEnd.out.normalized_files


TrinityNormalizeReads_DoubleEnd.out.normalized_files.view { meta, file ->
    "Normalized Double End File: Sample ID: ${meta.id}, File Name: $file.name | Path: $file"
}
}

//Take look this!!
//Normalized Double End File: Sample ID: all, File Name: [left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq, right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq] | Path: [/Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/work/e7/25404dede3fc0aa04868df4abc13f7/results_trinity/insilico_read_normalization_altogether/left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq, /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/work/e7/25404dede3fc0aa04868df4abc13f7/results_trinity/insilico_read_normalization_altogether/right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq]


//In terms of Single end
Channel.empty().set { ch_normalized_single_end_files }

if (params.single_end_sample) {

    ch_samples_single_end = ch_filtered_reads
        .filter { meta, path ->
            meta.single_end == true
        }

    //simultaneously run Trinity to expedite the run time
    TRINITY_NORMALIZATION_PARALLEL_SingleEnd(ch_samples_single_end )

    // collect() to collect all  transcript_fasta
    ch_samples_single_end = TRINITY_NORMALIZATION_PARALLEL_SingleEnd.out.transcript_fastq.collect()

    ch_versions = ch_versions.mix(TRINITY_NORMALIZATION_PARALLEL_SingleEnd.out.versions)


    ch_samples_single_end
          .flatMap()
            .map {
                 item ->
            if(item instanceof Map) {
            item['id'] = 'all_single'
            return item
            } else {
            return item
            }
            }
            .collate(2)
            .groupTuple()
            .map {
            meta, fastq ->
            [ meta, fastq.flatten() ]
            }
        .set { ch_inputfor_single_TrinityNormalization }

    // for single end
    ch_inputfor_single_TrinityNormalization.view { txt ->
        "Single End Sample Text Content: $txt"
    }

    TrinityNormalizeReads_SingleEnd(ch_inputfor_single_TrinityNormalization)

    ch_normalized_single_end_files = TrinityNormalizeReads_SingleEnd.out.normalized_files
    TrinityNormalizeReads_SingleEnd.out.normalized_files.view { meta, file ->
        "Normalized Single End File: Sample ID: ${meta.id}, File Name: $file.name | Path: $file"
    }
    ch_versions = ch_versions.mix(TrinityNormalizeReads_SingleEnd.out.versions)
}


// Using mix to combine the two channels to create ch_filtered_reads channel
ch_filtered_reads = ch_normalized_double_end_files.mix(ch_normalized_single_end_files)

//since we do not use trinity normalization by read set,we only can get one fastq.gz file.
if(params.trinity_normalization_by_read_set) {
    ch_filtered_reads = ch_filtered_reads.map { meta, file ->
        meta.single_end = true
        return [meta, file]
    }
}

// Giving ch_filtered_reads an alias to adapt to the input pattern of STAR
ch_filtered_reads_for_star = ch_filtered_reads

// View the contents of ch_filtered_reads_for_star
ch_filtered_reads_for_star.view { meta, file ->
    "Filtered Reads for STAR: Sample ID: ${meta.id}, Single_end: ${meta.single_end}, Standedness:  ${meta.strandedness}. File Name: $file.name | Path: $file"}

//Filtered Reads for STAR: Sample ID: all, Single_end: true, Standedness:  reverse. File Name: single.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq | Path: /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/work/a4/cd7d209b3906a76e5fea0b2387b77e/results_trinity/insilico_read_normalization_altogether/single.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq
//Filtered Reads for STAR: Sample ID: all, Single_end: false, Standedness:  reverse. File Name: [left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq, right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq] | Path: [/Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/work/0d/30707dbf1d9d65bbc1eb4c76c3af87/results_trinity/insilico_read_normalization_altogether/left.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq, /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/work/0d/30707dbf1d9d65bbc1eb4c76c3af87/results_trinity/insilico_read_normalization_altogether/right.norm.fq_ext_all_reads.normalized_K25_maxC200_minC1_maxCV10000.fq]
if (params.double_end_sample) {
ch_versions = ch_versions.mix(TrinityNormalizeReads_DoubleEnd.out.versions)
}
if (params.single_end_sample) {
ch_versions = ch_versions.mix(TrinityNormalizeReads_SingleEnd.out.versions)
}



if (params.fastqc_after_trinity)
{
        FASTQC_AFTER_TRINITY (ch_filtered_reads)
        fastqc_html = FASTQC_AFTER_TRINITY.out.html
        fastqc_zip  = FASTQC_AFTER_TRINITY.out.zip
        ch_versions = ch_versions.mix(FASTQC_AFTER_TRINITY.out.versions.first())} */





if (params.double_end_sample) {

ch_filtered_reads = ch_filtered_reads
    .filter { meta, path ->
        meta.single_end == false || meta.single_end == null  // null is for the case of undefined
    }
}


//In terms of Single end

if (params.single_end_sample) {

    ch_filtered_reads = ch_filtered_reads
        .filter { meta, path ->
            meta.single_end == true
        }

}

   // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index.map { [ [:], it ] },
            PREPARE_GENOME.out.gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        //please check subworkflow align_stat, you will found emit： bam            = BAM_SORT_STATS_SAMTOOLS.out.bam
        // and module/nf-core star/align/main.nf

        //we use ch_genome_bam        = ALIGN_STAR.out.bam_sorted
        // to replace original bam            = ALIGN_STAR.out.bam


        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final




        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

        //since '--outSAMtype BAM SortedByCoordinate' in module always have error 109
    // we use samtools to sort bam file
    // in config file  bam_csi_index  = false with_umi = false



    // SAMTOOLS_SORT_BAM(ch_genome_bam)

    // ch_genome_bam        = SAMTOOLS_SORT_BAM.out.bam
    // //ch_genome_bam_index = SAMTOOLS_SORT_BAM.out.cai

    // ch_versions = ch_versions.mix(SAMTOOLS_SORT_BAM.out.versions)






       // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
//  default with_umi = false
/*
        if (params.with_umi) {
            // Deduplicate genome BAM file before downstream analysis
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)


            // Co-ordinate sort, index and run stats on transcriptome BAM
            // !! important! only sort the transcriptome bam file
            BAM_SORT_STATS_SAMTOOLS (
                ch_transcriptome_bam,
                PREPARE_GENOME.out.fasta.map { [ [:], it ] }
            )
            ch_transcriptome_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
            ch_transcriptome_sorted_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

            // Deduplicate transcriptome BAM file before read counting with Salmon
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME (
                ch_transcriptome_sorted_bam.join(ch_transcriptome_sorted_bai, by: [0]),
                params.umitools_dedup_stats
            )



            // Name sort BAM before passing to Salmon
            SAMTOOLS_SORT (
                BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME.out.bam
            )

            // Only run prepare_for_rsem.py on paired-end BAM files
            SAMTOOLS_SORT
                .out
                .bam
                .branch {
                    meta, bam ->
                        single_end: meta.single_end
                            return [ meta, bam ]
                        paired_end: !meta.single_end
                            return [ meta, bam ]
                }
                .set { ch_umitools_dedup_bam }

            // Fix paired-end reads in name sorted BAM file
            // See: https://github.com/nf-core/rnaseq/issues/828
            UMITOOLS_PREPAREFORSALMON (
                ch_umitools_dedup_bam.paired_end
            )
            ch_versions = ch_versions.mix(UMITOOLS_PREPAREFORSALMON.out.versions.first())

            ch_umitools_dedup_bam
                .single_end
                .mix(UMITOOLS_PREPAREFORSALMON.out.bam)
                .set { ch_transcriptome_bam }
                ch_transcriptome_bam.view{ "Ready for bamsifter  Meta: ${it[0]}, Path: ${it[1]}" }
        }
 */
// BAMSIFTER normalizes the bam files in parallel.

    ch_bamstifer=BAMSIFTER(ch_genome_bam).normalized_bam.collect()

// put all the bamsifter outputs into a single channel
    ch_bamstifer
    .flatMap()
            .map {
                 item ->
            if(item instanceof Map) {
            item['id'] = 'all_double'
            return item
            } else {
            return item
            }
            }
            .collate(2)
            .groupTuple()
            .map {
            meta, fastq ->
            [ meta, fastq.flatten() ]
            }
            .set{ch_bamstifer_ready_samtools_merged}
            ch_bamstifer_ready_samtools_merged.view{ "Ready for Samtools_Merge  Meta: ${it[0]}, Path: ${it[1]}" }

        // Merge all the bam files using samtools merge
        //before merging, we need to sort the bam files
        //Forturnately, we have the bam files sorted in the previous step
        //and after we normalize the bam files, those normalized bam files are also sorted.


    SAMTOOLS_MERGE(ch_bamstifer_ready_samtools_merged,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fai.map { [ [:], it ] }
        )
     ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_MERGE.out.bam.view()

    //

    BAMSIFTER_NORMALIZATION_MERGED_BAM(SAMTOOLS_MERGE.out.bam)
    BAMSIFTER_NORMALIZATION_MERGED_BAM.out.normalized_bam.set{ch_genome_bam}
    ch_genome_bam.view{ "Ready for StringTie  Meta: ${it[0]}, Path: ${it[1]}" }






    //
    // MODULE: STRINGTIE_STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }

    //
    // MODULE: GFFCOMPARE
    //

    // combine fasta and fai into one channel，make sure format is tuple val(meta2), path(fasta), path(fai)
ch_combine_fasta_fai = PREPARE_GENOME.out.fasta.combine(PREPARE_GENOME.out.fai)
                          .map { fasta, fai -> [ [:], fasta, fai ] }
    // convert reference_gtf int a format  tuple val(meta3), path(reference_gtf)
ch_reference_gtf = PREPARE_GENOME.out.gtf.map { [ [:], it ] }

    GFFCOMPARE(STRINGTIE_STRINGTIE.out.transcript_gtf,ch_combine_fasta_fai,ch_reference_gtf)

    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())




    //we add this
    // MODULE: STRINGTIE_MERGE
    //The module output: tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    if (!params.skip_alignment && !params.skip_stringtie) {
        //modify the STRINGTIE_STRINGTIE output format to align the STRINGTIE_MERGE input format
        ch_stringtie_gtf_only = STRINGTIE_STRINGTIE.out.transcript_gtf.map { meta, gtf -> return gtf }

        STRINGTIE_MERGE (
            ch_stringtie_gtf_only,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions.first())
    }


    // update ch_transcriptome_bam, below code from prepare_genome, I just change  MAKE_TRANSCRIPTS_FASTA name and GTF_GENE_FILTER name

 /*     previous code like:
    else {
        ch_filter_gtf = GTF_GENE_FILTER ( ch_fasta, ch_gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    } */


    ch_filter_gtf = GTF_GENE_FILTER_FROM_NEW_GTF ( PREPARE_GENOME.out.fasta, STRINGTIE_MERGE.out.gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA_FROM_NEW_GTF ( PREPARE_GENOME.out.fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER_FROM_NEW_GTF.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA_FROM_NEW_GTF.out.versions)

      ch_new_salmon_index =  SALMON_INDEX_FROM_NEW_TRANSCRIPT_FASTA (

            PREPARE_GENOME.out.fasta,
             ch_transcript_fasta
        ).index

         ch_versions         = ch_versions.mix(SALMON_INDEX_FROM_NEW_TRANSCRIPT_FASTA.out.versions)














        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        // update previous PREPARE_GENOME.out.transcript_fasta, with new ch_transcript_fasta
        // update previous PREPARE_GENOME.out.gtf with new STRINGTIE_MERGE.out.gtf
        //when you use ch_transcript_fasta, will have error SAM file says target ENSDART00000152259 has length 934, but the FASTA file contains a sequence of length [1772 or 1771]
        // I calcualte  ENSDART00000152259  length in origin Danio_rerio.GRCz11.109.gtf file, the  length is 934
        /*

        Exon 1 Length=59147540−59147385+1=156
Exon 2 Length
=
59148361
−
59148231
+
1
=
131
Exon 2 Length=59148361−59148231+1=131
Exon 3 Length
=
59148500
−
59148457
+
1
=
44
Exon 3 Length=59148500−59148457+1=44
Exon 4 Length
=
59148808
−
59148759
+
1
=
50
Exon 4 Length=59148808−59148759+1=50
Exon 5 Length
=
59148970
−
59148887
+
1
=
84
Exon 5 Length=59148970−59148887+1=84
Exon 6 Length
=
59149173
−
59149052
+
1
=
122
Exon 6 Length=59149173−59149052+1=122
Exon 7 Length
=
59150220
−
59149874
+
1
=
347
Exon 7 Length=59150220−59149874+1=347
         */

         //when official pipeline use ch_dummy_file, the ch_dummy_file is empty, they just set PREPARE_GENOME.out.salmon_index  empty

        QUANTIFY_STAR_SALMON (
            // ch_transcriptome_bam,
            // ch_dummy_file,
            // PREPARE_GENOME.out.transcript_fasta,
            // STRINGTIE_MERGE.out.gtf,
            // true,
            // params.salmon_quant_libtype ?: ''


      // instead using bam as input we use raw fastq, since bam from old gtf.
      //we also need to update PREPARE_GENOME.out.salmon_index with ch_new_salmon_index  using ch_transcript_fasta

            ch_filtered_reads,
            ch_new_salmon_index,
            ch_dummy_file,
            STRINGTIE_MERGE.out.gtf,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_versions = ch_versions.mix(QUANTIFY_STAR_SALMON.out.versions)

        /* if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_STAR_SALMON (
                QUANTIFY_STAR_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_STAR_SALMON.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_STAR_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_STAR_SALMON.out.versions)
        } */
    }
    //this purple barcket is for star_align



    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with RSEM
    // in nextflow config file aligner = stat_salmon
  /*   ch_rsem_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_rsem') {
        QUANTIFY_RSEM (
            ch_filtered_reads,
            PREPARE_GENOME.out.rsem_index,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = QUANTIFY_RSEM.out.bam
        ch_genome_bam_index  = QUANTIFY_RSEM.out.bai
        ch_samtools_stats    = QUANTIFY_RSEM.out.stats
        ch_samtools_flagstat = QUANTIFY_RSEM.out.flagstat
        ch_samtools_idxstats = QUANTIFY_RSEM.out.idxstats
        ch_star_multiqc      = QUANTIFY_RSEM.out.logs
        ch_rsem_multiqc      = QUANTIFY_RSEM.out.stat
        if (params.bam_csi_index) {
            ch_genome_bam_index = QUANTIFY_RSEM.out.csi
        }
        ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_RSEM (
                QUANTIFY_RSEM.out.merged_counts_gene,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_aligner_pca_multiqc        = DESEQ2_QC_RSEM.out.pca_multiqc
            ch_aligner_clustering_multiqc = DESEQ2_QC_RSEM.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_RSEM.out.versions)
        }
    }


    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_filtered_reads,
            PREPARE_GENOME.out.hisat2_index.map { [ [:], it ] },
            PREPARE_GENOME.out.splicesites.map { [ [:], it ] },
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = FASTQ_ALIGN_HISAT2.out.summary
        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        // default  with_umi  = fasle

        if (params.with_umi) {
            BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                params.umitools_dedup_stats
            )
            ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)
        }


    }





    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_multiqc
            .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        ch_genome_bam
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam }

        ch_genome_bam_index
            .join(ch_percent_mapped, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_genome_bam_index }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_mapped_reads[meta.id] = true
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    pass_mapped_reads[meta.id] = false
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
    }

    //
    // MODULE: Run Preseq
    // skip_qc = false skip_preseq true

   /*  ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }
 */

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_markduplicates && !params.with_umi) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fai.map { [ [:], it ] }
        )
        ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
        if (params.bam_csi_index) {
            ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }






    //
    // MODULE: Feature biotype QC using featureCounts
    //
    ch_featurecounts_multiqc = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc && !params.skip_biotype_qc && biotype) {
 //update our  genome reference gtf file after applying STRINGTIE_MERGE module
    //using  STRINGTIE_MERGE.out.gtf to replace PREPARE_GENOME.out.gtf
        STRINGTIE_MERGE
            .out
            .gtf
            .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype

        //update our  genome reference gtf file after applying STRINGTIE_MERGE module
    //using  STRINGTIE_MERGE.out.gtf to replace PREPARE_GENOME.out.gtf
        ch_genome_bam
            .combine(STRINGTIE_MERGE.out.gtf)
            .combine(biotype_in_gtf)
            .filter { it[-1] }
            .map { it[0..<it.size()-1] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header_multiqc
        )
        ch_featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    //
    if (!params.skip_alignment && !params.skip_bigwig) {

        BEDTOOLS_GENOMECOV (
            ch_genome_bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            BEDTOOLS_GENOMECOV.out.bedgraph_forward,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD.out.versions)

        BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            BEDTOOLS_GENOMECOV.out.bedgraph_reverse,
            PREPARE_GENOME.out.chrom_sizes
        )
    }

    //
    // MODULE: Downstream QC steps
    //
    ch_qualimap_multiqc           = Channel.empty()
    ch_dupradar_multiqc           = Channel.empty()
    ch_bamstat_multiqc            = Channel.empty()
    ch_inferexperiment_multiqc    = Channel.empty()
    ch_innerdistance_multiqc      = Channel.empty()
    ch_junctionannotation_multiqc = Channel.empty()
    ch_junctionsaturation_multiqc = Channel.empty()
    ch_readdistribution_multiqc   = Channel.empty()
    ch_readduplication_multiqc    = Channel.empty()
    ch_fail_strand_multiqc        = Channel.empty()
    ch_tin_multiqc                = Channel.empty()
    if (!params.skip_alignment && !params.skip_qc) {
        if (!params.skip_qualimap) {
                  //update our  genome reference gtf file after applying STRINGTIE_MERGE module
    //using  STRINGTIE_MERGE.out.gtf to replace PREPARE_GENOME.out.gtf
            QUALIMAP_RNASEQ (
                ch_genome_bam,
                STRINGTIE_MERGE.out.gtf.map { [ [:], it ] }
            )
            ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results
            ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())
        }
             //update our  genome reference gtf file after applying STRINGTIE_MERGE module
    //using  STRINGTIE_MERGE.out.gtf to replace PREPARE_GENOME.out.gtf
        if (!params.skip_dupradar) {
            DUPRADAR (
                ch_genome_bam,
                STRINGTIE_MERGE.out.gtf
            )
            ch_dupradar_multiqc = DUPRADAR.out.multiqc
            ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
        }

        if (!params.skip_rseqc && rseqc_modules.size() > 0) {
            BAM_RSEQC (
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                PREPARE_GENOME.out.gene_bed,
                rseqc_modules
            )
            ch_bamstat_multiqc            = BAM_RSEQC.out.bamstat_txt
            ch_inferexperiment_multiqc    = BAM_RSEQC.out.inferexperiment_txt
            ch_innerdistance_multiqc      = BAM_RSEQC.out.innerdistance_freq
            ch_junctionannotation_multiqc = BAM_RSEQC.out.junctionannotation_log
            ch_junctionsaturation_multiqc = BAM_RSEQC.out.junctionsaturation_rscript
            ch_readdistribution_multiqc   = BAM_RSEQC.out.readdistribution_txt
            ch_readduplication_multiqc    = BAM_RSEQC.out.readduplication_pos_xls
            ch_tin_multiqc                = BAM_RSEQC.out.tin_txt
            ch_versions = ch_versions.mix(BAM_RSEQC.out.versions)

            ch_inferexperiment_multiqc
                .map {
                    meta, strand_log ->
                        def inferred_strand = WorkflowRnaseq.getInferexperimentStrandedness(strand_log, 30)
                        pass_strand_check[meta.id] = true
                        if (meta.strandedness != inferred_strand[0]) {
                            pass_strand_check[meta.id] = false
                            return [ "$meta.id\t$meta.strandedness\t${inferred_strand.join('\t')}" ]
                        }
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = [
                            "Sample",
                            "Provided strandedness",
                            "Inferred strandedness",
                            "Sense (%)",
                            "Antisense (%)",
                            "Undetermined (%)"
                        ]
                        WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_fail_strand_multiqc }
        }
    }

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon') {
                   //update our  genome reference gtf file after applying STRINGTIE_MERGE module
    //using  STRINGTIE_MERGE.out.gtf to replace PREPARE_GENOME.out.gtf
    //when official pipeline use ch_dummy_file, the ch_dummy_file is empty, they just set transcript_fasta empty

       /*  QUANTIFY_SALMON (
            ch_filtered_reads,
            PREPARE_GENOME.out.salmon_index,
            ch_dummy_file,
            STRINGTIE_MERGE.out.gtf,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_salmon_multiqc = QUANTIFY_SALMON.out.results
        ch_versions = ch_versions.mix(QUANTIFY_SALMON.out.versions)

        if (!params.skip_qc & !params.skip_deseq2_qc) {
            DESEQ2_QC_SALMON (
                QUANTIFY_SALMON.out.counts_gene_length_scaled,
                ch_pca_header_multiqc,
                ch_clustering_header_multiqc
            )
            ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_SALMON.out.pca_multiqc
            ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_SALMON.out.dists_multiqc
            ch_versions = ch_versions.mix(DESEQ2_QC_SALMON.out.versions)
        } */
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'),
            ch_multiqc_logo.collect().ifEmpty([]),
            ch_fail_trimming_multiqc.collectFile(name: 'fail_trimmed_samples_mqc.tsv').ifEmpty([]),
            //ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
            ch_fail_strand_multiqc.collectFile(name: 'fail_strand_check_mqc.tsv').ifEmpty([]),
            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),
            ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
            ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
            ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
            ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_aligner_pca_multiqc.collect().ifEmpty([]),
            ch_aligner_clustering_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
            ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
            ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
            ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
            ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
            ch_readduplication_multiqc.collect{it[1]}.ifEmpty([]),
            ch_tin_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()

    }

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
    }

    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS

To create the `StringTieMerge` branch on your remote server and synchronize it with your GitHub repository, follow these steps:

1. **Navigate to your repository directory**:
   ```bash
   cd /path/to/your/repo
   ```

2. **Ensure your git repository is linked to the remote**:
   ```
   git remote -v
   ```
   This should display the URL of your GitHub repository. If it doesn't, you'll need to add it:
   ```
   git remote add origin YOUR_GITHUB_REPO_URL

   use this to update your url
   git remote set-url origin https://github.com/dxu104/rnaseq_transcriptome_update.git
   ```

3. **Fetch all updates from GitHub for all branches**:
   ```
   git fetch origin
   ```

4. **Switch to the `StringTieMerge` branch**. If this branch does not exist locally, you'll need to create and switch to it:
   ```
   git checkout -b StringTieMerge origin/StringTieMerge
   git checkout -b  Bamsifter origin/Bamsifter
   git checkout -b  Bamsifter_Merge origin/Bamsifter_Merge
    git checkout -b update_transcript_fasta origin/update_transcript_fasta
    git checkout -b axolotltest origin/axolotltest
    git checkout -b axolotl_test_module_update origin/axolotl_test_module_update
   ```

By now, you should be on the `StringTieMerge` branch on your remote server, and it should be synchronized with the same branch in your GitHub repository.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// git merge origin/master or git merge origin/developbrach


// -stub-run
//cat nextflow.out
////tail nextflow.out
//ps aux |grep nextflow
// top  ctrl+S to lock the screen
//ps -u dxu
//kill -9  1234

// /compbio/referece/axolotl-omics/AmexT_v47.chr.unscaffolded.CherryGFP.gtf

//move the input file and json to the local laptop
//scp -r dxu@random.mdibl.org:/compbio/scratch/dxu/zfTest/ /Users/xudecheng/Library/Mobile\ Documents/com~apple~CloudDocs/MDIBL/RNAseq_TrinityNormalization/launch_dir/
//  scp -r dxu@random.mdibl.org:/compbio/scratch/dxu/newrnaseq/workdir/11/  /Users/dxu/whymerge_soslow/OnRamdon_output_workdir_cop_fromRandom
// scp  dxu@random.mdibl.org:/home/dxu/jfs_ssh.key  /Users/dxu/
//move the input file and json from local to the random
//aws s3 cp /compbio/analysis/JamesGodwin/jgodwin_001.nf_core_genome/star_salmon/stringtie/3A.coverage.gtf s3://biocore-data/JamesGodwin/jgodwin_001/gtfFile/
 //scp -r dxu@random.mdibl.org:/compbio/analysis/JamesGodwin/jgodwin_001.nf_core_genome/star_salmon/stringtie/3A.transcripts.gtf s3://biocore-data/JamesGodwin/jgodwin_001/gtfFile/
//scp -r /Users/dxu/MDI/RNAseq_TrinityNormalization/rnaseq/nextflow.AWSBatch.config dxu@random.mdibl.org:/compbio/scratch/dxu/newrnaseq/TestUnkownerroor

//copy s3 to random do not work we need download from s3 to local and then upload to random
// aws s3 cp s3://biocore-data/external/ensembl/release-109/danio_rerio.genome.fa /Users/dxu/MDI/RNAseq_TrinityNormalization/testMergetool
// scp -r /Users/dxu/MDI/RNAseq_TrinityNormalization/testMergetool/danio_rerio.genome.fa dxu@random.mdibl.org:/compbio/scratch/dxu/testMergetool

// scp -r /Users/dxu/MDI/RNAseq_TrinityNormalization/testMergetool dxu@random.mdibl.org:/compbio/scratch/dxu/testMergetool

//scp -r /Users/xudecheng/Library/Mobile\ Documents/com~apple~CloudDocs/MDIBL/RNAseq_TrinityNormalization/launch_dir/* dxu@random.mdibl.org:/compbio/scratch/dxu/newrnaseq/launch_dir/

//use -bg to run in the background https://www.nextflow.io/docs/latest/cli.html?highlight=bg


//copy file from random to aws s3 bucket

//aws s3 cp /compbio/reference/axolotl-omics/AmexG_v6_chr_unscaffolded.fa s3://biocore-data/JamesGodwin/jgodwin_001/axolotlGTFAndFasta/
//aws s3 cp /compbio/reference/axolotl-omics/AmexT_v47.chr.unscaffolded.CherryGFP.gtf s3://biocore-data/JamesGodwin/jgodwin_001/axolotlGTFAndFasta/



//Runner HeadNode
//copy file from HeadNode to local mac


//Copy file from headnode to S3
//This is original salmon_Index from REPARE_GENOME:SALMON_INDEX:  aws s3 sync /mnt/jfs/nextflow/dxu/axolotl2samples_Nov1st2023/9e/24ca52b68c1a08e40844f98ce7e3d4/salmon s3://biocore-data/JamesGodwin/jgodwin_001/AxolotlSalmonIndexFromREPARE_GENOME_SALMON_INDEX/

//  This is new salmon_Index after stringTie Merge:  aws s3 sync /mnt/jfs/nextflow/dxu/axolotl2samples_Nov1st2023/50/305d14ec4744f205ef1af91d49a1c7/salmon  s3://biocore-data/JamesGodwin/jgodwin_001/AxolotlSalmonIdex/

//  This is StarIndex:  aws s3 sync /mnt/jfs/nextflow/dxu/axolotl2samples_Nov1st2023/2c/2af926e34ca767e3b97a61747dfb39/star s3://biocore-data/JamesGodwin/jgodwin_001/AxolotlStarIndex/
// sync error aws s3 sync /mnt/jfs/nextflow/dxu/axolotl2samples_Nov20th2023_NoEoption_StringTie/50/893aef581150fb0833478268f65b13 s3://biocore-data/JamesGodwin/jgodwin_001/Error/StrandneitherPlus_Minues/
// scp -i /Users/dxu/jfs_admin_multiuser_ssh.key dxu@52.55.106.191:/mnt/jfs/nextflow/dxu/mmcloud.config.template /Users/dxu/MDI/RNAseq_TrinityNormalization/launch_dir/axolotl2samples/

// scp -i /Users/dxu/jfs_admin_multiuser_ssh.key dxu@52.55.106.191:/mnt/jfs/nextflow/dxu/axolotl2samples_Nov1st2023/ec/664b1b514bdc19e0114cbc64aac702 /Users/dxu/Documents/axolotStringTieMergeSuccess/outPutformFisrtSalmon_Index
// scp -i /Users/dxu/jfs_admin_multiuser_ssh.key dxu@52.55.106.191:/mnt/jfs/nextflow/dxu/axolotl2samples_Nov1st2023/35/1c49e8* /Users/dxu/Documents/axolotStringTieMergeSuccess/resultsFromGffcompare/
//scp -i /Users/dxu/jfs_admin_multiuser_ssh.key -r /Users/dxu/MDI/RNAseq_TrinityNormalization/launch_dir root@ec2-52-55-106-191.compute-1.amazonaws.com:/mnt/jfs/nextflow/dxu

// !!!this is the command to run the pipeline on AWSBATCH no memverge
// nextflow run main.nf -profile docker -c nextflow.AWSBatch.config -with-tower -work-dir s3://mdibl-nextflow-work/dxu/Bamsifter_why_merge_so_slow_AWSBatch_no_MemVerge -params-file ../launch_dir/zfTestAWSBatch/zf_params_AWSBatch.json
// nextflow run nf-core/fetchngs -profile docker test_full -c nextflow.AWSBatch.config --outdir s3://mdibl-dxu/TestWhyBatchNotWork -work-dir s3://mdibl-nextflow-work/dxu/TestWhyBatchNotWork  -with-tower
//
//Axolotl aws batch ondemand

//
// nextflow run main.nf -profile docker -c nextflow.AWSBatch.config -with-tower -work-dir s3://mdibl-nextflow-work/dxu/axolotl2samples_10-30-23_BatchIOPS -params-file ../launch_dir/axolotl2samples/parameter.json -resume

// axolotl mutliuse jfs memverge
//ssh -i /Users/dxu/jfs_admin_multiuser_ssh.key dxu@52.55.106.191
//  cd rnaseq_transcriptome_update
//tmux new -s your session name

//nextflow run main.nf -profile docker -c ../launch_dir/axolotl2samples/mmcloud.config  -params-file ../launch_dir/axolotl2samples/parameter.json -with-tower -resume
//OnDemand
//nextflow run main.nf -c ../launch_dir/axolotl2samples/mmcloud_OnDemand.config  -params-file ../launch_dir/axolotl2samples/parameter.json -with-tower -resume

//partialOnDemand
//nextflow run main.nf  -c ../launch_dir/axolotl2samples/mmcloud_PartialOnDemand.config  -params-file ../launch_dir/axolotl2samples/parameter.json -with-tower -resume



//on local mac
 //nextflow run main.nf -profile docker   -work-dir /Users/dxu/whymerge_soslow/localworkdir -params-file ../launch_dir/zfTestlocal/zf_paramslocal.json

//on random
//nextflow run main.nf -profile docker  -with-tower  -work-dir /compbio/scratch/dxu/newrnaseq/workdir -params-file ../launch_dir/zfTestRandom/zf_paramslocal.json

//memverge
//nextflow run main.nf -profile docker -c ../launch_dir/zfTestMemverge/float.config -with-tower  -params-file ../launch_dir/zfTestMemverge/zf_params_memverge.json  -resume -bg

//memverge ondeman
//nextflow run main.nf -profile docker -c ../launch_dir/zfTestMemvergeOndemand/float_ondemand.config -with-tower  -params-file ../launch_dir/zfTestMemvergeOndemand/zf_params_memvergeOndemand.json -resume

//memverge ondeman smallest test file
// nextflow run main.nf -profile test,docker -c ../launch_dir/smallestTestMemvergeOndemand/float_ondemand.config -with-tower  -params-file ../launch_dir/smallestTestMemvergeOndemand/zf_params_memvergeOndemand.json  -resume



//memverge ondeman two samples delete transmit without normalization by read set
//nextflow run main.nf -profile docker -c ../launch_dir/zfTestMemverge2samples/float.config -with-tower  -params-file ../launch_dir/zfTestMemverge2samples/zf_params_memverge.json -resume

//nextflow run main.nf -profile docker -c ../mmcloud.config  -params-file ../zf_params_memverge.json -with-tower -resume

//memverge ondeman two samples delete transmit and with normalization by read set
//nextflow run main.nf -profile docker -c ../launch_dir/zfTestMemverge2samplesnormalize_by_read_set/float.config -with-tower  -params-file ../launch_dir/zfTestMemverge2samplesnormalize_by_read_set/zf_params_memverge.json

//nohup nextflow run main.nf -profile test,docker -c nextflow.AWSBatch.config -with-tower --gene_prefix='AM-MIDBLv00003' -work-dir s3://mdibl-nextflow-work/dxu/zfish18files_09-18-23/ -params-file ../launch_dir/zfTest/zf_params.json --input ../launch_dir/zfTest/zfSamples.csv double_end_sample = true -resume &> nextflow.out&
//cat /proc/meminfo
// cd /compbio/scratch/dxu/newrnaseq/rnaseq
// cd /compbio/scratch/dxu/newrnaseq/rnaseq_copy
// cd /compbio/scratch/dxu/newrnaseq/rnaseq_smallest_test_delete_read_by_set
//tmux attach -t or tmux a -t
//tmux new -s your session name
//tmux rename-session -t preivousname new_name
//tmux kill-session -t your session name
//cat modules/local/TrinityNormalization/trinity_normalization.nf
//command to verify after bamsifter bam file is sorted or not.
//samtools view -H /Users/xudecheng/Downloads/all_double.reads1.bam | grep 'SO:'
//–outSAMtype: type of output. Default is BAM Unsorted; STAR outputs unsorted Aligned.out.bam file(s). “The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file cannot be directly used with downstream software such as HTseq, without the need of name sorting.” We therefore prefer the option BAM SortedByCoordinate
//sortmerna extremly slow even for 2mb test input
// Shift + option + A
//cd /compbio/scratch/dxu/testMergetool
//samtools merge  SL* -f -o all.bam --threads 15 --reference danio_rerio.genome.fa
//samtools view all.bam | parallel -j 32 --pipe wc -l | awk '{s+=$1} END {print s}'

// cat ~/.ssh/id_rsa.pub
// eval "$(ssh-agent -s)"
// ssh-add ~/.ssh/id_rsa
