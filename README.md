![CI](https://github.com/WEHIGenomicsRnD/qc-pipe/actions/workflows/main.yml/badge.svg)

# QC-pipe

A [snakemake](https://snakemake.readthedocs.io) pipeline for generating QC metrics from sequencing data.

### Installation ###

The only prerequisite is [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install snakemake, you will need to install a Conda-based Python3 distribution. For this, [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) is recommended. Once mamba is installed, snakemake can be installed like so:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Now activate the snakemake environment (you'll have to do this every time you want to run the pipeline):

```
conda activate snakemake
```

Now clone the repository:

```
git clone https://github.com/WEHIGenomicsRnD/qc-pipe.git
cd qc-pipe
```

### Testing ###

You can test the pipeline via:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 1 --directory .test
```

### Configuration ###

The configuration file is found under `config/config.yaml` and the config file for FastQ Screen is found under `config/fastq_screen.conf`. Please carefully go through these settings.

### Running ###

Place your fastq files in format of `{sample}_R[1|2].fastq.gz` under the the directory specific in your `config.yaml` file (`fastq` by default). Now run the pipeline as follows:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --cores 1
```

If you want to submit your jobs to the cluster using SLURM, use the following to run the pipeline:

```
conda activate snakemake
snakemake --use-conda --conda-frontend mamba --profile slurm --jobs 8 --cores 24
```

### Output ###

The pipeline will generate all results under a `results` directory. The most relevant directories are:

- `results/multiqc/multiqc_report.html` -- contains collated QC data in the form of an html report.
