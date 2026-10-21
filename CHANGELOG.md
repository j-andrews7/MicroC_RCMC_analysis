# j-andrews7/microcrcmc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

Nextflow port of the Hansen lab Micro-C/RCMC Snakemake pipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- nf-core template based Nextflow pipeline (nf-schema parameter/samplesheet validation, nf-core/configs institutional profiles such as `stjude`).
- Local modules reproducing every Snakemake rule with the original commands and tool versions; identical output folder layout and file names.
- Parallel alignment of FASTQ chunks (`--split_fastq`, nf-core `seqkit/split2`), re-joined in read order before sorting.
- Reuse of an existing bwa-mem2 index next to the FASTA, `--bwamem2_index` and `--save_reference` (nf-core `bwamem2/index`).
- Singularity/Docker containers (Seqera Community containers built from the original conda environment, plus a mustache container).
- Loop calling with mustache on condition-merged maps (`--mustache_resolutions`, per chromosome).
- Optional differential loop calling with diff_mustache (`--diff_loops`, `--contrasts`).
- FastQC and MultiQC (including pairtools deduplication statistics).
- `test`, `test_rcmc`, `test_diff` and `test_full` profiles and `bin/compare_outputs.py` to compare outputs with the original pipeline.

### `Changed`

- The Snakemake implementation moved to `legacy/snakemake/`.
- Samplesheet columns renamed to `sample,fastq_1,fastq_2,replicate,lane` (CSV).
- Logs are written to `<outdir>/logs/` instead of `./logs/`.
- RCMC balancing uses the job's CPUs instead of a hard-coded `--nproc 60` (results unchanged).

### `Dependencies`

| Tool      | Version |
| --------- | ------- |
| bwa-mem2  | 2.3     |
| pairtools | 1.1.3   |
| cooler    | 0.10.4  |
| bedops    | 2.4.42  |
| bedtools  | 2.31.1  |
| mustache  | 1.3.3   |
| seqkit    | 2.13.0  |
| fastqc    | 0.12.1  |
| multiqc   | 1.35    |