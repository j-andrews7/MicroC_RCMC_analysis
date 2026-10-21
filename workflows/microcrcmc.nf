/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome'
include { FASTQ_TO_PAIRS         } from '../subworkflows/local/fastq_to_pairs'
include { PAIRS_TO_MCOOLS        } from '../subworkflows/local/pairs_to_mcools'
include { MUSTACHE_LOOPS         } from '../subworkflows/local/mustache_loops'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_microcrcmc_pipeline'
include { parseResolutions       } from '../subworkflows/local/utils_nfcore_microcrcmc_pipeline'
include { getMustacheChromosomes } from '../subworkflows/local/utils_nfcore_microcrcmc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MICROCRCMC {

    take:
    ch_samplesheet // channel: [ meta, [ fastq_1, fastq_2 ] ] read in from --input
    ch_contrasts   // channel: [ condition1, condition2 ] for differential loop calling
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    def ch_chromsizes = channel.value(file(params.chromsizes, checkIfExists: true))
    def min_resolution = parseResolutions(params.resolutions).min()

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })
    }

    //
    // SUBWORKFLOW: bwa-mem2 index + RCMC blacklist
    //
    PREPARE_GENOME(
        params.fasta,
        params.bwamem2_index,
        params.chromsizes,
        params.regions
    )

    //
    // SUBWORKFLOW: align (chunked), parse, sort, merge lanes + dedup, merge replicates
    //
    FASTQ_TO_PAIRS(
        ch_samplesheet,
        PREPARE_GENOME.out.index,
        ch_chromsizes,
        params.assembly,
        params.mapq_min,
        params.split_fastq as Integer
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TO_PAIRS.out.dedup_stats.map { _meta, stats -> stats })

    //
    // SUBWORKFLOW: pairix + cooler cload + cooler zoomify/balance for replicates and conditions
    //
    PAIRS_TO_MCOOLS(
        FASTQ_TO_PAIRS.out.rep_pairs,
        FASTQ_TO_PAIRS.out.condition_pairs,
        ch_chromsizes,
        PREPARE_GENOME.out.blacklist,
        params.assembly,
        min_resolution,
        params.resolutions
    )

    //
    // SUBWORKFLOW: mustache loop calling (+ optional diff_mustache) on condition-merged maps
    //
    def call_loops = !params.skip_loops
    def diff_loops = params.diff_loops || params.contrasts
    if (call_loops || diff_loops) {
        MUSTACHE_LOOPS(
            PAIRS_TO_MCOOLS.out.condition_mcools,
            ch_contrasts,
            parseResolutions(params.mustache_resolutions),
            getMustacheChromosomes(params),
            call_loops,
            diff_loops as boolean
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'microcrcmc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'microcrcmc'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
