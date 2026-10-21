//
// Subworkflow with functionality specific to the j-andrews7/microcrcmc pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run j-andrews7/MicroC_RCMC_analysis -profile <singularity/docker/.../institute> --input samplesheet.csv --fasta genome.fa.gz --chromsizes genome.chrom.sizes --assembly <name> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    def samplesheet_rows = validateInputSamplesheet(
        samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    )
    def conditions = samplesheet_rows.collect { row -> row[0].condition }.unique()

    //
    // Validate pipeline-specific parameters and build the list of differential comparisons
    //
    validatePipelineParams(params)
    def contrast_list = getContrasts(params, conditions)

    ch_samplesheet = channel.fromList(samplesheet_rows)
    ch_contrasts   = channel.fromList(contrast_list)

    emit:
    samplesheet = ch_samplesheet
    contrasts   = ch_contrasts
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate the samplesheet rows and add the metadata needed to group lanes -> replicates -> conditions.
// Returns a list of [ meta, [ fastq_1, fastq_2 ] ] in samplesheet order.
//
def validateInputSamplesheet(rows) {
    def errors = []
    rows.each { row ->
        def (meta, fastq_1, fastq_2) = row
        if (fastq_1.name == fastq_2.name) {
            errors << "fastq_1 and fastq_2 must have different file names (${meta.condition}, ${meta.rep}, ${meta.lane}): ${fastq_1.name}"
        }
    }
    def seen = [] as Set
    rows.each { row ->
        def meta = row[0]
        def key = "${meta.condition}_${meta.rep}_${meta.lane}".toString()
        if (!seen.add(key)) {
            errors << "Duplicate sample/replicate/lane combination in samplesheet: ${meta.condition}, ${meta.rep}, ${meta.lane}"
        }
    }
    if (errors) {
        error("Please check input samplesheet ->\n  " + errors.join("\n  "))
    }

    def lanes_per_rep = rows.countBy { row -> "${row[0].condition}\t${row[0].rep}".toString() }
    def reps_per_condition = rows.collect { row -> [row[0].condition, row[0].rep] }
        .unique()
        .countBy { cond_rep -> cond_rep[0] }

    return rows.withIndex().collect { row, idx ->
        def (meta, fastq_1, fastq_2) = row
        def new_meta = meta + [
            id          : "${meta.condition}_${meta.rep}_${meta.lane}".toString(),
            single_end  : false,
            lane_order  : idx,
            n_lanes     : lanes_per_rep["${meta.condition}\t${meta.rep}".toString()],
            n_reps      : reps_per_condition[meta.condition],
            fastq1_name : fastq_1.name,
            fastq2_name : fastq_2.name,
        ]
        [new_meta, [fastq_1, fastq_2]]
    }
}

//
// Integer resolutions from a comma-separated string
//
def parseResolutions(String res_string) {
    return res_string.toString().tokenize(',').collect { res -> res.trim() as Integer }
}

//
// Validate pipeline-specific parameter combinations not expressible in the JSON schema
//
def validatePipelineParams(pipeline_params) {
    def errors = []
    if (!pipeline_params.fasta && !pipeline_params.bwamem2_index) {
        errors << "Either --fasta or --bwamem2_index must be provided."
    }
    def resolutions = parseResolutions(pipeline_params.resolutions)
    def min_res = resolutions.min()
    def bad_res = resolutions.findAll { res -> res % min_res != 0 }
    if (bad_res) {
        errors << "All --resolutions must be integer multiples of the smallest resolution (${min_res}); offending: ${bad_res.join(',')}"
    }
    if (!pipeline_params.skip_loops || pipeline_params.diff_loops || pipeline_params.contrasts) {
        def mustache_res = parseResolutions(pipeline_params.mustache_resolutions)
        def missing = mustache_res.findAll { res -> !(res in resolutions) }
        if (missing) {
            errors << "--mustache_resolutions must be a subset of --resolutions; not found: ${missing.join(',')}"
        }
        if (pipeline_params.mustache_chromosomes) {
            def known = file(pipeline_params.chromsizes).readLines().findAll { line -> line.trim() }.collect { line -> line.trim().tokenize()[0] }
            def unknown = pipeline_params.mustache_chromosomes.toString().tokenize(',').collect { chrom -> chrom.trim() }.findAll { chrom -> !(chrom in known) }
            if (unknown) {
                errors << "--mustache_chromosomes contains chromosomes not present in --chromsizes: ${unknown.join(',')}"
            }
        }
    }
    if (errors) {
        error("Invalid parameters ->\n  " + errors.join("\n  "))
    }
}

//
// Differential loop comparisons: from --contrasts if given, otherwise all pairwise combinations of conditions
// (in samplesheet order) when --diff_loops is set. Returns a list of [ condition1, condition2 ].
//
def getContrasts(pipeline_params, conditions) {
    if (!pipeline_params.diff_loops && !pipeline_params.contrasts) {
        return []
    }
    def contrasts = []
    if (pipeline_params.contrasts) {
        if (!pipeline_params.diff_loops) {
            log.info "--contrasts was provided, so differential loop calling (--diff_loops) is enabled."
        }
        contrasts = samplesheetToList(pipeline_params.contrasts, "${projectDir}/assets/schema_contrasts.json")
            .collect { row -> [row[0].toString(), row[1].toString()] }
        def errors = []
        contrasts.each { c1, c2 ->
            if (c1 == c2) {
                errors << "condition1 and condition2 must differ: ${c1}"
            }
            [c1, c2].each { cond ->
                if (!(cond in conditions)) {
                    errors << "Contrast condition '${cond}' is not a sample in the samplesheet (${conditions.join(', ')})"
                }
            }
        }
        if (errors) {
            error("Please check the contrasts file ->\n  " + errors.join("\n  "))
        }
    } else {
        if (conditions.size() < 2) {
            error("--diff_loops requires at least two conditions (samples) in the samplesheet; found: ${conditions.join(', ')}")
        }
        conditions.eachWithIndex { c1, i ->
            conditions.eachWithIndex { c2, j ->
                if (i < j) {
                    contrasts << [c1, c2]
                }
            }
        }
    }
    return contrasts
}

//
// Chromosomes to run mustache / diff_mustache on
//
def getMustacheChromosomes(pipeline_params) {
    if (pipeline_params.mustache_chromosomes) {
        return pipeline_params.mustache_chromosomes.toString().tokenize(',').collect { chrom -> chrom.trim() }
    }
    def exclude = (pipeline_params.mustache_exclude_chroms ?: '').toString().tokenize(',').collect { chrom -> chrom.trim() }
    def chroms = file(pipeline_params.chromsizes).readLines()
        .findAll { line -> line.trim() }
        .collect { line -> line.trim().tokenize() }
        .findAll { fields -> (fields[1] as Long) > (pipeline_params.mustache_min_chrom_size as Long) && !(fields[0] in exclude) }
        .collect { fields -> fields[0] }
    if (pipeline_params.regions) {
        def region_chroms = file(pipeline_params.regions).readLines()
            .findAll { line -> line.trim() && !line.startsWith('#') && !line.startsWith('track') && !line.startsWith('browser') }
            .collect { line -> line.trim().tokenize()[0] }
            .unique()
        chroms = chroms.findAll { chrom -> chrom in region_chroms }
    }
    if (!chroms) {
        error("No chromosomes left for mustache loop calling; check --mustache_chromosomes / --mustache_exclude_chroms / --regions.")
    }
    return chroms
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "bwa-mem2 (Vasimuddin et al. 2019),",
            "pairtools (Open2C et al. 2024),",
            "pairix (Lee et al. 2022),",
            "cooler (Abdennur & Mirny 2020),",
            params.regions ? "BEDOPS (Neph et al. 2012), BEDtools (Quinlan & Hall 2010)," : "",
            (!params.skip_loops || params.diff_loops || params.contrasts) ? "mustache (Roayaei Ardakany et al. 2020)," : "",
            params.split_fastq ? "SeqKit (Shen et al. 2016)," : "",
            params.skip_fastqc ? "" : "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Vasimuddin Md, Misra S, Li H, Aluru S (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. IEEE IPDPS, 314-324. doi: 10.1109/IPDPS.2019.00041</li>",
            "<li>Open2C, Abdennur N, Fudenberg G, Flyamer IM, Galitsyna AA, Goloborodko A, Imakaev M, Venev SV (2024). Pairtools: from sequencing data to chromosome contacts. PLoS Comput Biol 20(5): e1012164. doi: 10.1371/journal.pcbi.1012164</li>",
            "<li>Lee S, Bakker CR, Vitzthum C, Alver BH, Park PJ (2022). Pairs and Pairix: a file format and a tool for efficient storage and retrieval for Hi-C read pairs. Bioinformatics 38(6): 1729-1731. doi: 10.1093/bioinformatics/btab870</li>",
            "<li>Abdennur N, Mirny LA (2020). Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics 36(1): 311-316. doi: 10.1093/bioinformatics/btz540</li>",
            params.regions ? "<li>Neph S, Kuehn MS, Reynolds AP, et al. (2012). BEDOPS: high-performance genomic feature operations. Bioinformatics 28(14): 1919-1920. doi: 10.1093/bioinformatics/bts277</li>" : "",
            params.regions ? "<li>Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6): 841-842. doi: 10.1093/bioinformatics/btq033</li>" : "",
            (!params.skip_loops || params.diff_loops || params.contrasts) ? "<li>Roayaei Ardakany A, Gezer HT, Lonardi S, Ay F (2020). Mustache: multi-scale detection of chromatin loops from Hi-C and Micro-C maps using scale-space representation. Genome Biol 21: 256. doi: 10.1186/s13059-020-02167-0</li>" : "",
            params.split_fastq ? "<li>Shen W, Le S, Li Y, Hu F (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. doi: 10.1371/journal.pone.0163962</li>" : "",
            params.skip_fastqc ? "" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".").replaceAll(" +", " ")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
