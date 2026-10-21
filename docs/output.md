# j-andrews7/microcrcmc: Output

All outputs are written to `--outdir`. The first five folders are identical in name and content to the outputs of the original Snakemake pipeline.

```
<outdir>/
├── beds/             # RCMC only: chroms.bed, blacklist.bed
├── dedup/            # per replicate: deduplicated pairs, unmapped/duplicate pairs, stats, pairix index
├── repmerge/         # per condition: replicate-merged deduplicated pairs + pairix index
├── repcools/         # per replicate: .cool (base resolution) and balanced multi-resolution .mcool
├── conditioncools/   # per condition: .cool and balanced .mcool
├── loops/            # mustache loop calls (new)
├── diffloops/        # diff_mustache differential loops (new, optional)
├── fastqc/           # FastQC reports (new)
├── multiqc/          # MultiQC report (new)
├── logs/             # stderr of every step (the original pipeline wrote these to ./logs)
├── genome/           # bwa-mem2 index (only with --save_reference)
└── pipeline_info/    # Nextflow execution reports, software versions, parameters
```

## Pairs and deduplication

`dedup/` (per replicate `<sample>_<replicate>`):

- `<sample>_<replicate>_nodups.pairs.gz`: deduplicated pairs of all lanes of the replicate (bgzipped, sorted, [pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md), columns `readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2` plus walk columns from `pairtools parse2 --expand`; read IDs are dropped).
- `<sample>_<replicate>_nodups.pairs.gz.px2`: pairix index.
- `<sample>_<replicate>_unmapped.pairs.gz`: pairs with at least one unmapped or multi-mapped (MAPQ < `--mapq_min`) side.
- `<sample>_<replicate>_dups.pairs.gz`: PCR/optical duplicates.
- `<sample>_<replicate>.dedup.stats`: `pairtools` statistics (shown in the MultiQC report).

`repmerge/` (per condition):

- `<sample>_merged_nodups.pairs.gz` (+ `.px2`): deduplicated pairs of all replicates of the condition, merged without further deduplication.

## Contact maps

`repcools/` and `conditioncools/`:

- `<name>.cool`: contact map at the smallest resolution of `--resolutions`.
- `<name>.mcool`: multi-resolution contact map (`--resolutions`), balanced at every resolution (`weight` column). In RCMC mode, bins outside the capture regions are blacklisted from balancing (NaN weights).

where `<name>` is `<sample>_<replicate>` (replicates) or `<sample>_merged` (conditions).

`beds/` (RCMC only):

- `chroms.bed`: whole-chromosome intervals from `--chromsizes`.
- `blacklist.bed`: everything outside the capture regions (`--regions`), used for balancing.

### Viewing maps in HiGlass

`.mcool` files can be viewed in [HiGlass](https://higlass.io/) (Kerpedjiev et al., Genome Biology 19:125, 2018). For example, with a local HiGlass Docker container:

```bash
mkdir -p higlass/data higlass/tmp
docker pull higlass/higlass-docker
docker run --detach --publish XXXX:80 \
    --volume /full/path/to/higlass/data:/data \
    --volume /full/path/to/higlass/tmp:/tmp \
    --name higlass_container_name higlass/higlass-docker
cp results/conditioncools/WT_merged.mcool higlass/tmp/
docker exec higlass_container_name python higlass-server/manage.py ingest_tileset \
    --filename /tmp/WT_merged.mcool --filetype cooler --datatype matrix --project-name my_mcools
```

Then open `http://<hostname>:XXXX/app` and add the map from the `+` menu. Loop calls from `loops/` can be converted to BEDPE (columns 1-6) and loaded as 2D annotations.

## Loops

`loops/`:

- `<sample>_merged_<resolution>bp.mustache.tsv`: mustache loop calls on the condition-merged map at each `--mustache_resolutions` resolution, concatenated over chromosomes. Columns: `BIN1_CHR BIN1_START BIN1_END BIN2_CHROMOSOME BIN2_START BIN2_END FDR DETECTION_SCALE`.

`diffloops/` (with `--diff_loops` / `--contrasts`), per comparison `<condition1>_vs_<condition2>` and resolution:

- `<c1>_vs_<c2>_<res>bp.loop1`: loops found in `<condition1>`.
- `<c1>_vs_<c2>_<res>bp.loop2`: loops found in `<condition2>`.
- `<c1>_vs_<c2>_<res>bp.diffloop1`: loops present in `<condition1>` that are weakened or lost in `<condition2>`.
- `<c1>_vs_<c2>_<res>bp.diffloop2`: loops present in `<condition2>` that are weakened or lost in `<condition1>`.

All have the same columns as the mustache output.

## QC

- `fastqc/`: FastQC reports per lane (`<sample>_<replicate>_<lane>_{1,2}_fastqc.{html,zip}`).
- `multiqc/multiqc_report.html`: FastQC, pairtools deduplication statistics (pair types, cis/trans, duplicates, distance distributions), software versions and the parameters used.

## Logs

`logs/` contains the stderr of each step, mirroring the original pipeline's `logs/` folder:

- `alignment/<sample>_<replicate>_<lane>[.part_NNN].log`: bwa-mem2 + pairtools parse2 (per chunk); `<...>.sort.log`: per-lane sort.
- `dedup/`, `merge/`, `repcools/`, `conditioncools/`, `repmcools/`, `conditionmcools/`: the corresponding steps.
- `mustache/`, `diff_mustache/`: per condition/comparison, resolution and chromosome.

## Pipeline information

`pipeline_info/` contains Nextflow's execution report, timeline, trace and DAG, the parameters used (`params_<timestamp>.json`) and the versions of all tools (`microcrcmc_software_mqc_versions.yml`).
