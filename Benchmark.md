# Benchmarking TKSM
Welcome to the benchmarking branch of TKSM!
This document will help you reproduce all the results presented in the manuscript.

- [Benchmarking TKSM](#benchmarking-tksm)
  - [Setup](#setup)
    - [Environment](#environment)
    - [Data](#data)
  - [Running](#running)
    - [Regular run](#regular-run)
    - [Piped](#piped)
  - [Results](#results)


## Setup
The main dependency that you will need to have installed is [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) or, much prefereblally, [Mamba](https://github.com/conda-forge/miniforge#mambaforge).
If you choose to install Conda instead of Mamba, replace `mamba` in all the commands with `conda`.
Also, expect a much slower installation.


### Environment
Clone the repository and its submodules:
```bash
git clone https://github.com/vpc-ccg/tksm
cd tksm
git checkout paper
git submodule update --init --recursive
```

Create a new Conda environment satisfying all the dependencies:
```bash
mamba env create -f Benchmark_env.yaml
mamba activate tksm_bench
```

All the following commands assume that you are in the `tksm` directory and that the `tksm_bench` environment is activated.

Install TKSM:
```bash
make -j8
./install.sh
```

### Data

Create the data directories:
```bash
mkdir -p data/refs
mkdir -p data/samples
```

Download the samples:
```bash
wget http://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/fastq/SGNex_MCF7_directRNA_replicate2_run2/SGNex_MCF7_directRNA_replicate2_run2.fastq.gz -O data/samples/MCF7-sgnex.fastq.gz

wget https://figshare.com/ndownloader/files/40862291?private_link=f66abe6156696e14ced6 -O data/samples/N1.fastq.gz
```

Download the references:
```bash
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -O data/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz -O data/refs/Homo_sapiens.GRCh38.108.chr.gtf.gz

wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O data/refs/Homo_sapiens.GRCh38.cdna.chr.fa.gz

wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz -O data/refs/genomicSuperDups.txt.gz
```

Note: `wget` can be slow especially with the references.
Consider using [Axel](https://github.com/axel-download-accelerator/axel) which you can install using [Conda](https://anaconda.org/conda-forge/axel). 

Unzip the references:
```bash
gunzip data/refs/genomicSuperDups.txt.gz
gunzip data/refs/Homo_sapiens.GRCh38.108.chr.gtf.gz
gunzip data/refs/Homo_sapiens.GRCh38.cdna.chr.fa.gz
gunzip data/refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Download the 10x Genomics cellular barcode whitelist:
```bash
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz -O data/refs/3M-february-2018.txt.gz
```

## Running
The benchmark pipelines are specified by the following files:

- `Snakefile`: TKSM's Snakefile specifying TKSM module rules and tracks its memory and time using GNU Time `v1.9`. 
- `config.yaml`: TKSM's config file. This is not used by this benchmarking pipelines.
- `Benchmark.smk`: Imports all the rules from `Snakefile` and adds [Trans-Nanosim](https://github.com/bcgsc/NanoSim) rules. It also has all the rules for generating the result plots.
- `Benchmark_config.yaml`: Specifies the expirement pipelines, input/output for all the pipelines. Edit this to change where things are output or to add new experiments!
- `config_piped.yaml`: A special version of the `Benchmark_config.yaml` that allows for running only TKSM pipelines in piped mode. The config is run by `Snakefile`, and NOT by `Benchmark.smk`.


### Regular run
To run the benchmarking pipelines, run the following command:
```bash
snakemake -s Benchmark.smk -j <threads>
```

This will run all the pipelines specified in `Benchmark_config.yaml` and output the results in `output` directory.

### Piped
First make sure that `time.tsv` file as been created:
```bash
snakemake --configfile config_piped.yaml -U make_time
```

To benchmark the different pipelines, each experiment pipeline must be run separately.
Therefore, open `config_piped.yaml` and uncomment the experiment you want to run and make sure that all the other experiments are commented out.
By default the `P_TKSM_bulk` experiment is the only uncommented experiment.

Once you have uncommented the experiment you want to run, run the following command:
```bash
d=$(date +"%Y-%m-%d %H:%M:%S"); \
$(which time) -f "$d,P_TKSM_bulk,,,%e,%U,%M,,%S,True" -o >(awk 'BEGIN{FS=","} {printf "%s,%s,%s,%s,%.2f,%.2f,%.2f,%d,%.2f,%s\n",$1,$2,$3,$4,$5/60,$6/60,$7/(1024*1024),$8,$9/60,$10}' >> output/time.tsv) \
snakemake --configfile config_piped.yaml -p -j <threads>
```

This first captures the current time and stores it in the variable `d`.
Then it runs the GNU `time` command (not the bash keyword!) on the whole Snakemake pipeline.
The `-f` flag specifies the format of the output.
Please edit the `-f` value to match the experiment you are running (in this case its `P_TKSM_bulk`).
The benchmarked GNU Time values are then piped to `awk` which parses the output and appends it to `output/time.tsv` file after formatting them (e.g. converting memory use from `KB` to `GB`).
Note that if you changed the output path in the `config_piped.yaml` file, you will need to change the `awk` command's output file to match the new output path.

Repeat the above steps for each experiment you want to run.

## Results
The results are all under the `output` directory (unless you changed the output path in the `Benchmakr_config.yaml` or `config_piped.yaml` files).
The main files of interest are:

- `time.tsv`: Contains the time and memory usage of each pipeline (headers are in the file).
- `plots/`: Contains all the plots generated by `Benchmark.smk` (e.g. `plots/substitution/MCF7-sgnex.TKSM_bulk.nanosim_bulk.pdf`).

Note that the plot files start with a prefix composed of the `"."` concatenation of the input sample name (e.g. `MCF7-sgnex`) followed by the experiments names (e.g. `TKSM_bulk.nanosim_bulk`).

The preprocessing, intermediate, and output files of TKSM and Trans-Nanosim are all included in the `output` directory.