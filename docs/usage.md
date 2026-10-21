# j-andrews7/microcrcmc: Usage

> Parameters can be given on the command line (`--param value`) or, preferably, in a YAML/JSON file passed with `-params-file params.yml`. Run `nextflow run j-andrews7/MicroC_RCMC_analysis --help` (or `--help_full`) to list all parameters.

## Introduction

The pipeline processes paired-end Micro-C or RCMC FASTQ files into deduplicated pairs, multi-resolution balanced contact maps (`.mcool`) and chromatin loop calls. It replicates the original Snakemake pipeline (kept in [`legacy/snakemake/`](../legacy/snakemake/)) step by step and adds parallelism, containers, QC and loop calling.

## Samplesheet input

Use `--input` to pass a comma-separated (`.csv`) or tab-separated (`.tsv`) file with a header row and the following columns:

| Column | Description |
| --- | --- |
| `sample` | Condition name (e.g. `WT`, `KO`). Replicates of the same `sample` are merged into `<sample>_merged`. No spaces or `/`. |
| `fastq_1` | Path to the gzipped FASTQ file of read 1 (`.fastq.gz` / `.fq.gz`). |
| `fastq_2` | Path to the gzipped FASTQ file of read 2. Must have a different file name than `fastq_1`. |
| `replicate` | Replicate label (any text, e.g. `rep1`). Required even if there is only one replicate. |
| `lane` | Lane/run label (any text, e.g. `lane1`). Required even if there is only one lane. |

Each `sample` + `replicate` + `lane` combination must be unique. Example: [`assets/samplesheet.csv`](../assets/samplesheet.csv).

```csv
sample,fastq_1,fastq_2,replicate,lane
WT,/data/WT_rep1_L001_R1.fastq.gz,/data/WT_rep1_L001_R2.fastq.gz,rep1,lane1
WT,/data/WT_rep1_L002_R1.fastq.gz,/data/WT_rep1_L002_R2.fastq.gz,rep1,lane2
WT,/data/WT_rep2_L001_R1.fastq.gz,/data/WT_rep2_L001_R2.fastq.gz,rep2,lane1
KO,/data/KO_rep1_L001_R1.fastq.gz,/data/KO_rep1_L001_R2.fastq.gz,rep1,lane1
```

Output names follow the original pipeline: `<sample>_<replicate>` for replicate-level files and `<sample>_merged` for condition-level files.

Relative paths in the samplesheet are resolved relative to the directory Nextflow is launched from; absolute paths are safest.

## Reference files

| Parameter | Description |
| --- | --- |
| `--fasta` | Genome FASTA (may be gzipped). Used to build the bwa-mem2 index if none is provided. If the five bwa-mem2 index files (`<fasta>.0123`, `.amb`, `.ann`, `.bwt.2bit.64`, `.pac`) already exist next to the FASTA, as created by the original Snakemake pipeline, they are used directly. |
| `--bwamem2_index` | Directory containing exactly one pre-built bwa-mem2 index (e.g. `<outdir>/genome/bwamem2` from a previous run with `--save_reference`). |
| `--save_reference` | Save the index built by the pipeline to `<outdir>/genome/`. Building a human index takes several hours and ~90 GB of memory, so save it once and reuse it. |
| `--chromsizes` | Two-column, tab-separated chromosome sizes file (chromosome, length) with only the major chromosomes, naturally sorted (chr1, chr2, ..., chr10, ..., chrX, chrY, chrM). See [`demo/hg19.sorted.chrom.sizes`](../demo/hg19.sorted.chrom.sizes). |
| `--assembly` | Assembly name stored in the pairs/cool files (e.g. `hg38`). |

## Processing options

| Parameter | Default | Description |
| --- | --- | --- |
| `--resolutions` | `10000000,...,1000` | Comma-separated resolutions for `cooler zoomify`. All must be integer multiples of the smallest, which is used for `cooler cload`. |
| `--regions` | none | RCMC mode: BED file (chrom, start, end, no header) with the captured regions. Balancing then ignores everything outside these regions (blacklist written to `beds/`). |
| `--mapq_min` | `30` | Minimum MAPQ for `pairtools parse2`. |
| `--split_fastq` | `20000000` | Read pairs per alignment chunk; `0` aligns each lane as a single job. |
| `--skip_fastqc` | `false` | Skip FastQC. |
| `--skip_multiqc` | `false` | Skip MultiQC. |

### Alignment chunking and reproducibility

With `--split_fastq N`, each lane is split into chunks of `N` read pairs (`seqkit split2`), which are aligned and parsed in parallel. The parsed chunks are re-joined in their original order and sorted with one `pairtools sort` per lane. Because `pairtools sort` is stable and `bwa-mem2 mem -SP` does not use per-batch insert-size statistics, the resulting per-lane pairs are equivalent to aligning the whole lane in one job.

The only expected difference: when two alignments of a read have exactly the same score, bwa-mem2 breaks the tie with a hash of the read's index in the input file. Chunking changes those indices, so for a very small fraction of chimeric reads a different segment may be reported as primary. That can change the (capped) MAPQ of the supplementary segment, or which segments of a multi-ligation walk are paired. Set `--split_fastq 0` to reproduce the original pipeline exactly.

Measured on the demo data (3M read pairs aligned in 6 chunks vs. a single job, compared with `bin/compare_outputs.py`):

- 3 of 1,851,356 deduplicated pairs differed (0.0002%), e.g. a MAPQ of 33 instead of 60 on one side.
- A single contact moved between two pixels of the contact maps; totals and pixel counts were identical.
- Balancing weights changed by a median relative difference of ~1e-6 (Pearson r = 1.000000).

Pick a chunk size that gives alignment jobs of a convenient length on your cluster. Smaller chunks give more parallelism, but each job has to load the full bwa-mem2 index (on the order of a minute or two for a human genome).

## Loop calling (mustache)

Loops are called by default with [mustache](https://github.com/ay-lab/mustache) on every condition-merged `.mcool` (`<sample>_merged.mcool`), using the balancing weights. Each resolution x chromosome combination runs as a separate job, and the results are concatenated per condition and resolution in chromosome order.

| Parameter | Default | Description |
| --- | --- | --- |
| `--skip_loops` | `false` | Skip loop calling. |
| `--mustache_resolutions` | `1000,5000,10000` | Resolutions (bp) to call loops at; must be in `--resolutions`. |
| `--mustache_pt` | `0.1` | `-pt`, FDR threshold. |
| `--mustache_st` | `0.8` | `-st`, sparsity threshold. Lower values (e.g. `0.7`) are recommended for sparse data or very high resolutions. |
| `--mustache_sz` | mustache default (1.6) | `-sz`, sigma0. |
| `--mustache_octaves` | mustache default (2) | `-oc`. |
| `--mustache_iterations` | mustache default (10) | `-i`. |
| `--mustache_distance` | mustache default | `-d`, maximum loop distance (bp). mustache uses a resolution-dependent default capped at 2 Mb. |
| `--mustache_chromosomes` | automatic | Comma-separated chromosomes to use, overriding the automatic selection. |
| `--mustache_exclude_chroms` | `chrM,MT,M,chrMT` | Chromosomes never used. |

Automatic chromosome selection: all chromosomes in `--chromsizes` longer than 1 Mb (mustache's own whole-genome behaviour) that are not excluded. In **RCMC mode**, only chromosomes containing at least one capture region are used, because balancing weights are undefined (NaN) everywhere else.

Chromosomes without contacts produce an empty (header-only) loop file. If mustache fails on a chromosome for any other reason (non-signal exit), the pipeline logs a warning in `logs/mustache/` and writes an empty file for it instead of stopping (hidden parameter `--mustache_allow_failure false` disables this). Out-of-memory kills are always retried with more memory.

## Differential loops (diff_mustache)

Enable with `--diff_loops`. By default all pairwise comparisons between conditions are run (in samplesheet order: `A_vs_B`, `A_vs_C`, `B_vs_C`, ...). To choose specific comparisons, pass `--contrasts contrasts.csv` (this also enables differential loop calling):

```csv
condition1,condition2
WT,KO
WT,HET
```

For each comparison and each `--mustache_resolutions` value, `diff_mustache` runs per chromosome on `<condition1>_merged.mcool` and `<condition2>_merged.mcool` with the mustache parameters above plus `--diff_mustache_pt2` (default `0.1`, the FDR threshold for differential loops).

## Running on St. Jude HPCF

Use `-profile stjude` (loaded automatically from [nf-core/configs](https://github.com/nf-core/configs/blob/master/conf/stjude.config)). It submits every process to LSF (`short`/`standard`/`large_mem` queues chosen by time and memory), runs it in Singularity, caps jobs at 32 CPUs and 1 TB of memory, and retries jobs killed for exceeding memory. Containers are downloaded once into Nextflow's Singularity cache; set `NXF_SINGULARITY_CACHEDIR` to a shared location so later runs reuse them.

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache
NXF_VER=25.10.4 nextflow run j-andrews7/MicroC_RCMC_analysis -r main -profile stjude -params-file params.yml
```

## Resources

Default requests per job are defined in [`conf/base.config`](../conf/base.config) and are capped by the profile's `resourceLimits`:

| Process | CPUs | Memory | Time |
| --- | --- | --- | --- |
| `BWAMEM2_PARSE2` (per chunk) | 16 | 48 GB | 24 h |
| `PAIRTOOLS_SORT_LANE` (per lane) | 16 | 64 GB | 24 h |
| `PAIRTOOLS_MERGE_DEDUP` (per replicate) | 8 | 32 GB | 48 h |
| `PAIRTOOLS_MERGE_REPS` (per condition) | 8 | 16 GB | 24 h |
| `COOLER_CLOAD_*` | 4 | 32 GB | 48 h |
| `COOLER_ZOOMIFY_*` | 16 | 64 GB | 48 h |
| `MUSTACHE` (per condition x resolution x chromosome) | 8 | 32 GB | 12 h |
| `DIFF_MUSTACHE` (per comparison x resolution x chromosome) | 8 | 64 GB | 24 h |

Memory (and time) is multiplied by the attempt number on retry. Override with a custom config passed via `-c`:

```groovy
process {
    withName: 'BWAMEM2_PARSE2' {
        cpus   = 32
        memory = 64.GB
    }
}
```

Tool arguments can be extended through `ext.args` (see [`conf/modules.config`](../conf/modules.config)), e.g. `ext.args2` of `BWAMEM2_PARSE2` is appended to `pairtools parse2`.

## Comparing with the original pipeline

1. Run the legacy pipeline (see [`legacy/snakemake/README.md`](../legacy/snakemake/README.md)) and this pipeline on the same inputs and reference. Use `--split_fastq 0` to check exact equivalence, or the default to measure the (tiny) effect of chunking.
2. Compare the output folders with [`bin/compare_outputs.py`](../bin/compare_outputs.py), e.g. inside the pipeline's tools container, which has `numpy` and `h5py`:

   ```bash
   singularity exec \
       https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6ed82fac71a9112629fa056825a5f04819b2ec1b702f1d3c3ce90e3b1559097d/data \
       python bin/compare_outputs.py legacy/snakemake/outputs results
   ```

   The script checks that every legacy output exists, then compares pairs records (ignoring headers), dedup statistics, BED files, and the chromosomes, bins, pixels and balancing weights of every resolution of every `.cool`/`.mcool`. It exits with status 1 if anything differs. Use `--skip-pairs` to skip the (slow) pairs comparison for large datasets.

The `test_full` profile reproduces the original demo (demo reads, full hg19), so a legacy demo run can be compared directly with:

```bash
nextflow run j-andrews7/MicroC_RCMC_analysis -profile test_full,stjude --outdir nf_demo --fasta /path/to/indexed/hg19.fa.gz
```

## Core Nextflow arguments

- `-profile`: configuration profile(s), e.g. `stjude`, `singularity`, `docker`, `apptainer`, `conda`, `test`, `test_rcmc`, `test_diff`, `test_full`. Profiles can be combined (`-profile test,singularity`).
- `-resume`: reuse cached results of a previous run (the St. Jude profile deletes the work directory after a successful run).
- `-c`: additional config file, e.g. for resource overrides.
- `-params-file`: YAML/JSON file with pipeline parameters.
- `-r`: pipeline revision (branch, tag or commit) to run.
