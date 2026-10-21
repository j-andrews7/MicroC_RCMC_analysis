# MicroC_RCMC_analysis (legacy Snakemake pipeline)

> [!NOTE]
> This folder contains the **original Snakemake implementation** of the pipeline, kept for reference and for comparing results with the Nextflow port at the root of this repository (see the [main README](../../README.md)). The `Snakefile` is unchanged; only the demo paths in `config/config.yaml` and `samples.tsv` were updated (`../../demo/...`) so the workflow runs from this folder. Run all commands below from `legacy/snakemake/`, and download the demo genome into the repository's `demo/` folder. Outputs are written to `legacy/snakemake/outputs/`.

A processing pipeline for alignment, processing, and generation of contact maps for Micro-C and RCMC data.

## Setup
### Conda setup (can be skipped if you already have conda installed)
The pipeline uses conda environments to run all of the necessary tools for processing the data. If you do not have conda installed, install conda and mamba as follows (or by any other method appropriate for your system), answering `yes` to all asked questions:
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```
Restart your terminal and conda will be activated.

### Environment setup
In order to run snakemake, a conda environment containing snakemake is needed. First, ensure the base environment is activated:
```
conda activate base
```
Next, create the required environment (conda can also be used here if mamba is not installed):
```
mamba env create --name snakemake-env -c conda-forge -c bioconda python=3.12 snakemake=9.12
```
or
```
conda create --name snakemake-env -c conda-forge -c bioconda python=3.12 snakemake=9.12
```
Once the environment is made, activate it:
```
conda activate snakemake-env
```
### Preparing files to run
Create and move to a directory for your analysis:
```
mkdir /path/to/your/analysis/directory
cd /path/to/your/analysis/directory
```
Download the repository to the directory and move to the directory:
```
git clone https://github.com/j-andrews7/MicroC_RCMC_analysis
cd MicroC_RCMC_analysis/legacy/snakemake
```
---
---
## Running demo analysis
### Setup
The repository includes files for running a quick demo analysis on heavily downsampled data. The fastq files for this data contain reads mapped to the human chromosome 16 from the hg19 build (note that these only contain 3,000,000 reads total, so the resulting maps will be sparse!). All of the files included are already configured to run on the demo data, but a genome file is needed. This can be downloaded with the following command, or you can use your own file (see below):
```
wget --directory-prefix ../../demo https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
```
---
### Running demo analysis - Micro-C

Run the pipeline, setting an appropriate number of cores for your system - this number serves as the upper limit for what your processes will use. If setting this lower than 8, also change the value set for the `threads` parameter in `config/config.yaml` using a text editor. On the first run, this will create the conda environment needed to run the pipeline, which may take some time, and then index the genome file that was downloaded - this may take ~2 hours. Once you already have run both of these once, they should not need to be run again. If you already have a `bwa-mem2` indexed genome to use, you can use this instead by modifying the `genome` parameter in `config/config.yaml` to have the full path to the location of the genome file and its associated 
```
snakemake --cores 8 --sdm conda
```
---
### Running demo analysis - RCMC
RCMC analysis can be run in much the same way as Micro-C analysis. To run in RCMC mode, comment out the blank `regions` entry in `config/config.yaml` and uncomment `regions: ../../demo/capture_locus.bed`. This provides a `.bed` format file which indicates the position of a captured region (no capture was performed on the demo data, but the same procedure can be performed to run the pipeline this way). Once `config/config.yaml` has been modified, you can again run the demo analysis to produce a new `.mcool` file which only contains reads in the indicated region:
```
snakemake --cores 8 --sdm conda
```
---
### Visualisation of demo data
Once execution completes, you can upload your completed `.mcool` file (found in `outputs/conditioncools`) to a [HiGlass](https://higlass.io/) instance[^1] using the following steps:
#### Create HiGlass directories
These can be created anywhere - we'll create them inside the directory where snakemake ran to keep everything in the same place.
```
mkdir higlass
mkdir higlass/data
mkdir higlass/tmp
```
#### Download the latest version of Docker HiGlass and create a HiGlass container
Replace XXXX with a port on your server that you can access and change the container name to something suitable. Ensure that the paths to your data and tmp directories are full paths, not relative.
```
docker pull higlass/higlass-docker
docker run --detach \
           --publish XXXX:80 \
           --volume /path/to/your/analysis/directory/MicroC_RCMC_analysis/higlass/data:/data \
           --volume /path/to/your/analysis/directory/MicroC_RCMC_analysis/higlass/tmp:/tmp \
           --name higlass_container_name \
           higlass/higlass-docker
```
#### Upload data
This copies your `.mcool` into the upload directory for HiGlass and then uploads it, storing it in a folder in the HiGlass UI with the name given for `--project-name`.
```
cp outputs/conditioncools/WT_merged.mcool higlass/tmp/
docker exec higlass_container_name python higlass-server/manage.py ingest_tileset --filename /tmp/WT_merged.mcool --filetype cooler --datatype matrix --project-name demo_mcools
```
#### Access your data on HiGlass
Go to [hostname.org:XXXX/app](http://hostname.org:XXXX/app) to load your HiGlass instance, where hostname.org is replaced with the address of your server (e.g. example.mit.edu) or its IP address, and XXXX is replaced with the port number specified above. You can load your `.mcool` for visualisation from the '+' menu in the upper right corner in the centre square. To situate yourself in the genome, you can also add chromosome position markers (under chromosomes) and gene annotations (under Gene Annotations) by choosing the appropriate options for your genome build.

---
---
## Running your own analysis
To run analysis of your own data, you will need to prepare the following things:
1. A new `samples.tsv` file (can be named whatever you want). Follow the template of the file already provided for including your `.fastq` files. Lane and replicate numbers must be filled, even if you only have one (but the values can be any arbitrary text).
2. If the data is from an RCMC experiment, a `.bed` file containing the locations of each of your captured regions. Follow the formatting of the provided `capture_locus.bed` file.
3. The appropriate genome for your analysis, and corresponding chromosome sizes file (can be downloaded from UCSC, as with the genome). The chromosome sizes file must follow the format of the provided template `hg19.sorted.chrom.sizes` file - namely, only complete chromosomes, and all autosomes sorted in numerical order, followed by sex chromosomes and mitotic DNA.
4. An updated `config.yaml` file, in a `config/` directory (relative to the location of the `Snakefile`) containing the paths to the above mentioned files in the correct locations, and with the `assembly`, `resolutions`, `threads` and `mapqmin` fields updated as appropriate. If analysing RCMC, provide the path to a `.bed` file in the `regions` field, otherwise leave it blank.

Once all of these are prepared, you can run your analysis as before using:
```
snakemake --cores 8 --sdm conda
```

[^1]: Kerpedjiev et al. HiGlass: Web-based visual comparison and exploration of genome interaction maps. Genome Biology, 19:125 (2018)
