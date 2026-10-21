import uuid
import pandas as pd
import os

# Name of config file
configfile: "config/config.yaml"

samples = (
	pd.read_csv(config["samples"], sep="\t", dtype={"condition": str})
	.set_index("condition", drop=False)
	.sort_index()
)

# Need minimum bin size (as a string) to bin cools to
reslist = config["resolutions"].split(",")
minres = str(min([int(res) for res in reslist]))

# Determine whether balancing should be normal (genome-wide) or only for regions
if config["regions"] is None:
	balancearg = ""
else:
	balancearg = f"--balance-args '--nproc 60 --blacklist {config['outdir']}beds/blacklist.bed' "

def sort_tmp():
	# Make randomly generated folder name for temporary sorting
	sortdir_id = str(uuid.uuid4())
	return sortdir_id

rule all:
	input:
		expand(config["outdir"]+"repcools/{condition}_{rep}.mcool", zip, condition = samples.condition, rep = samples.rep),
		expand(config["outdir"]+"conditioncools/{condition}_merged.mcool", zip, condition = samples.condition)


def make_chrombed(chromsizes, outdir):
	chromsizesdf = pd.read_csv(chromsizes, sep = "\t", names = ["chrom", "end"], header = None)
	chromsizesdf.insert(1, "start", [0] * len(chromsizesdf))
	chromsizesdf.to_csv(outdir+"beds/chroms.bed", sep = "\t", header = False, index = False)

# Make blacklist for balancing (RCMC only) - first make chromsizes into a bed file
rule make_chrombed:
	input:
		chromsizes = config["chromsizes"]
	output:
		config["outdir"]+"beds/chroms.bed"
	run:
		make_chrombed(input.chromsizes, config["outdir"])

# Then make blacklist bed
rule make_blacklist:
	input:
		chroms = config["outdir"]+"beds/chroms.bed",
		regions = branch(
			config["regions"],
			then = config["regions"],
			otherwise = [],
			)
	output:
		config["outdir"]+"beds/blacklist.bed"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"bedops --partition <(bedtools sort -i {input.chroms}) "
		"<(bedtools sort -i {input.regions}) | bedops "
		"--not-element-of 1 - <(bedtools sort -i {input.regions}) > {output}"

# Index genome file if necessary
rule index_genome:
	input:
		config["genome"]
	conda:
		"env/rcmc_conda.yaml"
	output:
		multiext(config["genome"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
	shell:
		"bwa-mem2 index {input}"

rule make_sort_dir:
	# input:
	# 	outdir = config["outdir"]
	output:
		temp(directory(config["outdir"]+"{condition}_sorttemp"))
	shell:
		"mkdir -p {output}"

# Function to get fastq names
def get_fastqs(wildcards):
	condition_name = f"{wildcards.condition}"
	rep_num = f"{wildcards.rep}"
	lane_num = f"{wildcards.lane}"
	return {
		"fastq1": samples["fastq1"][(samples["condition"] == condition_name) & (samples["rep"] == rep_num) & (samples["lane"]  == lane_num)],
		"fastq2": samples["fastq2"][(samples["condition"] == condition_name) & (samples["rep"] == rep_num) & (samples["lane"]  == lane_num)],
	}

# Align fastqs and generate pairs files. Does not save bams to save space.
rule make_pairs:
	input:
		unpack(get_fastqs),
		genome = config["genome"],
		genome_index = config["genome"]+".amb",
		sorttemp = config["outdir"]+"{condition}_sorttemp"
	output:
		temp(config["outdir"]+"pairs/{condition}_{rep}_{lane}.pairs")
	threads:
		config["threads"]
	params:
		assembly = config["assembly"],
		chromsizes = config["chromsizes"],
		mapqmin = config["mapqmin"]
	log:
		"logs/alignment/{condition}_{rep}_{lane}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		# Alignment command
		"(bwa-mem2 mem -t {threads} -SP {input.genome} {input.fastq1} {input.fastq2} | pairtools "
		"parse2 --add-columns mapq -c {params.chromsizes} --assembly {params.assembly} "
		"--min-mapq {params.mapqmin} --drop-sam --drop-readid --expand --nproc-in "
		"{threads} | pairtools sort --tmpdir {input.sorttemp} --nproc {threads} -o "
		"{output} | cat) 2> {log}"

# def lane_count(wildcards):
# 	lane_list = list(samples[(samples["condition"] == wildcards.condition) & (samples["rep"] == wildcards.rep)]["lane"])
# 	return len(lane_list)

def get_pairs_for_lane_merge(wildcards):
	all_lanes = list(samples[
		(samples["condition"] == wildcards.condition) &
		(samples["rep"] == wildcards.rep)
	]["lane"])
	return expand(rules.make_pairs.output, 
	condition=wildcards.condition, 
	rep=wildcards.rep, 
	lane=all_lanes)

# Merge lanes for the same sample and remove duplicates
rule merge_and_dedup_lanes:
	input:
		pairs = get_pairs_for_lane_merge,
		sorttemp = config["outdir"]+"{condition}_sorttemp"
	output:
		dedup = config["outdir"]+"dedup/{condition}_{rep}_nodups.pairs.gz",
		unmapped = config["outdir"]+"dedup/{condition}_{rep}_unmapped.pairs.gz",
		dups = config["outdir"]+"dedup/{condition}_{rep}_dups.pairs.gz",
		stats = config["outdir"]+"dedup/{condition}_{rep}.dedup.stats"
	threads:
		config["threads"]
	log:
		"logs/dedup/{condition}_{rep}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(pairtools merge --tmpdir {input.sorttemp} --nproc {threads} {input.pairs} | pairtools "
		"dedup --max-mismatch 1 --mark-dups --output {output.dedup} --output-unmapped "
		"{output.unmapped} --output-dups {output.dups} --output-stats {output.stats} | cat) "
		"2> {log}"

# # Remove duplicates without merging, for samples with only one lane, or where lanes were merged in advance
# rule dedup:
# 	input:
# 		# config["outdir"]+"{condition}_{rep}.pairs"
# 		pairs = expand(rules.make_pairs.output, zip, condition = samples.condition, rep = samples.rep, lane = samples.lane)
# 	output:
# 		dedup = config["outdir"]+"{condition}_{rep}_nodups.pairs.gz",
# 		unmapped = config["outdir"]+"{condition}_{rep}_unmapped.pairs.gz",
# 		dups = config["outdir"]+"{condition}_{rep}_dups.pairs.gz",
# 		stats = config["outdir"]+"{condition}_{rep}.dedup.stats"
# 	threads:
# 		config["threads"]
# 	log:
# 		"logs/dedup/{condition}_{rep}.log"
# 	conda:
# 		"env/rcmc_conda.yaml"
# 	shell:
# 		"(pairtools dedup --max-mismatch 1 --mark-dups --output {output.dedup} "
# 		"--output-unmapped {unmapped.dedup} --output-dups {output.dups} "
# 		"--output-stats {output.stats} {input}) 2> {log}"

def get_pairs_for_rep_merge(wildcards):
	all_reps = list(set(list(samples[
		(samples["condition"] == wildcards.condition)
	]["rep"])))
	return expand(config["outdir"]+"dedup/{condition}_{rep}_nodups.pairs.gz", 
	condition = wildcards.condition, 
	rep = all_reps)

# Merge replicates for the same condition without removing duplicates
rule merge_reps_by_condition:
	input:
		# dedup = config["outdir"]+"{{condition}}_{rep}_nodups.pairs.gz",
		dedup = get_pairs_for_rep_merge,
		sorttemp = config["outdir"]+"{condition}_sorttemp"
	output:
		config["outdir"]+"repmerge/{condition}_merged_nodups.pairs.gz"
	threads:
		config["threads"]
	log:
		"logs/merge/{condition}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(pairtools merge --tmpdir {input.sorttemp} --nproc {threads} "
		"--output {output} "
		"{input.dedup}) 2> {log}"

# Index pairs files for each replicate
rule index_reps:
	input:
		rules.merge_and_dedup_lanes.output.dedup
	output:
		config["outdir"]+"dedup/{condition}_{rep}_nodups.pairs.gz.px2"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"pairix {input}"

# Can probably combine this rule with the previous and process all at once
# Also index merged replicates
rule index_conditions:
	input:
		rules.merge_reps_by_condition.output
	output:
		config["outdir"]+"repmerge/{condition}_merged_nodups.pairs.gz.px2"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"pairix {input}"

# Make cool files for replicates
rule rep_cools:
	input:
		pairs = rules.merge_and_dedup_lanes.output.dedup,
		index = rules.index_reps.output
	output:
		config["outdir"]+"repcools/{condition}_{rep}.cool"
	threads:
		config["threads"]
	params:
		assembly = config["assembly"],
		chromsizes = config["chromsizes"],
		minimumresolution = minres
	log:
		"logs/repcools/{condition}_{rep}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(bgzip -cd -@ {threads} {input.pairs} | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "
		"--assembly {params.assembly} {params.chromsizes}:{params.minimumresolution} - {output}) 2> {log}"

# Make cool files for conditions
rule condition_cools:
	input:
		pairs = rules.merge_reps_by_condition.output,
		index = rules.index_conditions.output
	output:
		config["outdir"]+"conditioncools/{condition}_merged.cool"
	threads:
		config["threads"]
	params:
		assembly = config["assembly"],
		chromsizes = config["chromsizes"],
		minimumresolution = minres
	log:
		"logs/conditioncools/{condition}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(bgzip -cd -@ {threads} {input.pairs} | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "
		"--assembly {params.assembly} {params.chromsizes}:{params.minimumresolution} - {output}) 2> {log}"

# Make mcools for replicates
rule rep_mcools:
	input:
		cool = rules.rep_cools.output,
		# blacklist = config["outdir"]+"beds/blacklist.bed" if config["regions"] else []
		blacklist = branch(
			config["regions"],
			then = config["outdir"]+"beds/blacklist.bed",
			otherwise = [],
			)
	output:
		config["outdir"]+"repcools/{condition}_{rep}.mcool"
	threads:
		config["threads"]
	params:
		resolutions = config["resolutions"],
		balance = balancearg
	log:
		"logs/repmcools/{condition}_{rep}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(cooler zoomify --nproc {threads} --balance {params.balance}--out {output} --resolutions {params.resolutions} {input.cool}) 2> {log}"

# Make mcools for conditions
rule condition_mcools:
	input:
		cool = rules.condition_cools.output,
		# blacklist = config["outdir"]+"beds/blacklist.bed" if config["regions"] else []
		blacklist = branch(
			config["regions"],
			then = config["outdir"]+"beds/blacklist.bed",
			otherwise = [],
			)
	output:
		config["outdir"]+"conditioncools/{condition}_merged.mcool"
	threads:
		config["threads"]
	params:
		resolutions = config["resolutions"],
		balance = balancearg
	log:
		"logs/conditionmcools/{condition}.log"
	conda:
		"env/rcmc_conda.yaml"
	shell:
		"(cooler zoomify --nproc {threads} --balance {params.balance}--out {output} --resolutions {params.resolutions} {input.cool}) 2> {log}"
