# MicroC_RCMC_analysis (Nextflow)

> [!WARNING]
> **This Nextflow port was written with an AI coding assistant, and its output has not been verified against the original pipeline.** It has only been run on the small demo dataset (hg19 chr16, 3M read pairs). No full-size dataset has been processed with it, and its results have not been compared with the original Snakemake pipeline's outputs on real data. Until that comparison is done (see [Reproducing the original results](#reproducing-the-original-results)), treat all outputs, including pairs, contact maps and loop calls, as unvalidated. The original pipeline is kept unchanged in [`legacy/snakemake/`](legacy/snakemake/).

A parallelised [Nextflow](https://www.nextflow.io/) pipeline for alignment, processing, contact-map generation and loop calling of **Micro-C** and **RCMC** (Region Capture Micro-C) data, built with [nf-core](https://nf-co.re/docs/nf-core-tools) tooling (nf-core template, nf-schema parameter/samplesheet validation, nf-core modules and institutional configs).

This is a port of the original [Hansen lab Snakemake pipeline](https://github.com/ahansenlab/MicroC_RCMC_analysis), which is kept unchanged in [`legacy/snakemake/`](legacy/snakemake/). The core outputs (pairs, deduplication statistics, `.cool`/`.mcool` contact maps, RCMC blacklists) are produced with **the same commands and tool versions** and are written to **the same file names and folders**. The port adds:

- **Parallel execution on HPC**: every lane, replicate and condition runs as its own job, and each lane can be split into chunks that are aligned in parallel (`--split_fastq`).
- **Singularity/Docker containers** for every step (no conda environment needed).
- **Loop calling** with [mustache](https://github.com/ay-lab/mustache) on the condition-merged maps (on by default), parallelised per chromosome and resolution.
- **Optional differential loop calling** with `diff_mustache` (`--diff_loops` / `--contrasts`).
- FastQC and a MultiQC report (including pairtools deduplication statistics).

## Pipeline overview

1. Split each lane's FASTQs into chunks (`seqkit split2`; `--split_fastq` read pairs per chunk).
2. Align + parse each chunk: `bwa-mem2 mem -SP | pairtools parse2 --add-columns mapq --min-mapq 30 --drop-sam --drop-readid --expand`.
3. Re-join chunks in their original read order and sort per lane (`pairtools sort`).
4. Merge lanes per replicate and deduplicate (`pairtools merge | pairtools dedup --max-mismatch 1 --mark-dups`).
5. Merge replicates per condition (`pairtools merge`).
6. Index pairs (`pairix`), build contact maps at the smallest resolution (`cooler cload pairs`), then zoomify and balance (`cooler zoomify --balance`). RCMC: balancing is restricted to the captured regions via a blacklist (`bedops`/`bedtools`).
7. Call loops on the condition-merged `.mcool` files with mustache (per resolution x chromosome), and optionally differential loops between conditions with diff_mustache.
8. FastQC on raw reads and a MultiQC report.

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (**25.10.4 or a later 25.10.x release is recommended**, see [Nextflow version](#nextflow-version)) and make Singularity/Apptainer (or Docker) available.
2. Prepare a samplesheet (CSV):

   ```csv
   sample,fastq_1,fastq_2,replicate,lane
   WT,/data/WT_rep1_L001_R1.fastq.gz,/data/WT_rep1_L001_R2.fastq.gz,rep1,lane1
   WT,/data/WT_rep1_L002_R1.fastq.gz,/data/WT_rep1_L002_R2.fastq.gz,rep1,lane2
   WT,/data/WT_rep2_L001_R1.fastq.gz,/data/WT_rep2_L001_R2.fastq.gz,rep2,lane1
   KO,/data/KO_rep1_L001_R1.fastq.gz,/data/KO_rep1_L001_R2.fastq.gz,rep1,lane1
   ```

   `sample` is the experimental condition (the `condition` column of the old `samples.tsv`); `replicate` and `lane` can be any labels. Rows with the same `sample` + `replicate` are merged across lanes, and replicates of the same `sample` are merged into a condition-level map.

3. Run:

   ```bash
   nextflow run j-andrews7/MicroC_RCMC_analysis \
       -profile singularity \
       --input samplesheet.csv \
       --fasta /path/to/hg38.fa.gz \
       --chromsizes /path/to/hg38.sorted.chrom.sizes \
       --assembly hg38 \
       --outdir results
   ```

   For RCMC data add `--regions capture_regions.bed`.

### St. Jude HPCF

The [St. Jude institutional profile](https://github.com/nf-core/configs/blob/master/conf/stjude.config) is loaded automatically from nf-core/configs. It submits jobs to LSF and runs them in Singularity, so no other profile is needed:

```bash
NXF_VER=25.10.4 nextflow run j-andrews7/MicroC_RCMC_analysis -profile stjude -params-file params.yml
```

with e.g. `params.yml`:

```yaml
input: samplesheet.csv
outdir: results
fasta: /path/to/hg38.fa.gz
chromsizes: /path/to/hg38.sorted.chrom.sizes
assembly: hg38
# regions: capture_regions.bed      # RCMC only
# diff_loops: true                  # optional differential loops (all pairwise)
# contrasts: contrasts.csv          # or: specific comparisons (condition1,condition2)
```

Run the Nextflow head job somewhere it can stay alive for the whole run (e.g. an interactive session or a long-running `bsub` job). The St. Jude profile sets `cleanup = true`, so work directories are deleted after a successful run and `-resume` only helps for failed/interrupted runs.

### Nextflow version

The pipeline requires Nextflow `>=25.10.4`. With Nextflow **26.10**, values given on the command line reach the pipeline as strings, so typed options such as `--skip_loops` or `--split_fastq 0` fail parameter validation. Either use Nextflow 25.10.x (e.g. `NXF_VER=25.10.4 nextflow run ...`) or pass options through a `-params-file` (YAML/JSON), which keeps their types.

## Moving from the Snakemake pipeline

| Snakemake `config.yaml` | Nextflow parameter | Notes |
| --- | --- | --- |
| `samples` (TSV: condition, fastq1, fastq2, rep, lane) | `--input` (CSV: sample, fastq_1, fastq_2, replicate, lane) | same information, renamed columns |
| `outdir` | `--outdir` | same sub-folder layout |
| `genome` | `--fasta` | an existing bwa-mem2 index next to the FASTA is reused; or use `--bwamem2_index` |
| `assembly` | `--assembly` | |
| `chromsizes` | `--chromsizes` | |
| `resolutions` | `--resolutions` | same default |
| `regions` | `--regions` | RCMC mode when set |
| `mapqmin` | `--mapq_min` | default 30 |
| `threads` | per-process resources in [`conf/base.config`](conf/base.config) | capped by your profile, e.g. 32 CPUs on St. Jude HPCF |

To convert an existing `samples.tsv`:

```bash
awk 'BEGIN { FS = "\t"; OFS = "," } NR == 1 { print "sample,fastq_1,fastq_2,replicate,lane"; next } { print $1, $2, $3, $4, $5 }' samples.tsv > samplesheet.csv
```

## Reproducing the original results

- Every core step runs the same command as the Snakemake rule it replaces, in a container built from the original conda environment (`bwa-mem2 2.3`, `pairtools 1.1.3`, `cooler 0.10.4`, `bedops 2.4.42`, `bedtools 2.31`).
- `bwa-mem2 mem -SP` output does not depend on the number of threads or the batch size, so alignment can use as many CPUs as are available.
- **FASTQ chunking (`--split_fastq`, default 20,000,000 read pairs)**: chunks are re-joined in their original order before the (stable) `pairtools sort`, so the per-lane pairs are equivalent to aligning the lane in one job. The one exception is bwa-mem2's tie-breaking between alignments with exactly equal scores, which is seeded by each read's index in the input file. For a very small fraction of chimeric reads this can change which segment is reported as primary (and therefore a capped MAPQ, or the pairing within a multi-ligation walk). On the demo data, 6 chunks vs. 1 job changed 3 of 1.85M deduplicated pairs (0.0002%) and moved a single contact in the contact maps. Use `--split_fastq 0` to reproduce the original pipeline exactly.
- Header lines of pairs files (which embed command lines, paths and thread counts) and the creation date stored in cool files always differ between runs, including between two runs of the original pipeline.
- RCMC balancing originally used a hard-coded `--nproc 60`; it now uses the CPUs allocated to the job, which does not change the balancing weights.

[`bin/compare_outputs.py`](bin/compare_outputs.py) compares a legacy output folder with a Nextflow output folder (pairs records, dedup statistics, beds, and bins/pixels/weights at every cool/mcool resolution). See [docs/usage.md](docs/usage.md#comparing-with-the-original-pipeline).

## Documentation

- [docs/usage.md](docs/usage.md): inputs, parameters, profiles, loop calling, differential loops, resources, testing.
- [docs/output.md](docs/output.md): description of all outputs and how to view them in HiGlass.

## Test profiles

| Profile | Description |
| --- | --- |
| `test` | Demo reads (3M chr16 read pairs) aligned to hg19 chr16 only, in chunks of 500,000 read pairs; loops on chr16 |
| `test_rcmc` | As `test`, in RCMC mode with the demo capture regions |
| `test_diff` | Two conditions (multiple lanes/replicates) built from the demo reads, with `diff_mustache` |
| `test_full` | The original Snakemake demo: demo reads, full hg19 reference (building the index needs ~90 GB of memory) |

```bash
nextflow run j-andrews7/MicroC_RCMC_analysis -profile test,singularity --outdir test_results
```

## Credits

The original Snakemake pipeline was written by the [Hansen lab](https://github.com/ahansenlab). The Nextflow port is maintained by Jared Andrews ([@j-andrews7](https://github.com/j-andrews7)).

## Citations

See [CITATIONS.md](CITATIONS.md) for the tools used by the pipeline. Please also cite mustache if you use the loop calls:

> Roayaei Ardakany A, Gezer HT, Lonardi S, Ay F. Mustache: multi-scale detection of chromatin loops from Hi-C and Micro-C maps using scale-space representation. Genome Biol 21, 256 (2020). https://doi.org/10.1186/s13059-020-02167-0

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x.
