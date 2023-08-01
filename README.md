This repository is intended to be used as a Snakemake workflow template.
It provides all the necessary files and inline documentation to assist with generating your own customise workflow.

# Running the Workflow

To run a Snakemake workflow, you will first need to have gone through the [Snakemake One-Time Setup] to install Snakemake in a conda environment.

```bash
SNAKEMAKE_VERSION='7.8.1'

# Activate the conda environment to make Snakemake
# available to you on the command line
conda activate \
  "snakemake_${SNAKEMAKE_VERSION}"

# Run the workflow using the sahmri-hpc profile
# in order to have jobs submitted to Slurm
export LOCAL_SCRATCH=/cancer/storage/SAGC/scratch
snakemake \
  --profile profiles/sahmri-hpc
```

# Pipeline Setup

Create a project-specific `config.yaml` file from the supplied template:

```bash
cp config_template.yaml config.yaml
```

## Reference sequence

Prepare your reference sequence by ensuring each sequence record occupies 2 lines and is bgzip compressed, ready for SAMtools indexing.
This can be done by passing it through `fasta_formatter` of FASTX-Toolkit and using `bgzip` of HTSlib.
The file should be named `*.fasta.gz` and reside in the `references/` directory.
e.g.:

```bash
# Prep reference genome
dos2unix < my.fasta \
  | fasta_formatter \
  | bgzip -@ 12 \
> my.fasta.gz

# Index reference genome
samtools faidx my.fasta.gz
```

Next, you should update the `REFERENCE` variable in `config.yaml` to match the reference genome file basename.

## Gene Models

Gene models are used in the STAR genome indexing step.
Prepare these by creating a GTF file with the same filename as the reference genome, except for using `.gtf` instead of `.gz` file extension.
A GTF file may need to be created from a GFF3 file using `gffread` from Cufflinks.

```bash
# Correctly sort a GFF3 file
gff3sort.pl --precise references/my.gff3 \
  | bgzip \
  > references/my.gff3.gz

# Index the GFF3 file
tabix --csi references/my.gff3.gz

# Create GTF file from GFF3 file
pigz -dcp2 references/my.gff3.gz \
  | gffread -T - -o references/my.gtf

# Create a symlink to ensure the GTF file has the same prefix as the reference genome FASTA file
ln -s my.gtf references/my.fasta.gtf
```

## Sequence Data

The Snakemake pipeline expects the sequence read data to be structured in a specific way.
This is so the pipeline can automatically detect and process files without the need for editing configuration files.

Firstly, QC'd reads are placed under the `qc_reads/` directory with the following directory structure:

```bash
qc_reads/
└── transcriptomic
    └── PE
        └── dataset
            ├── sample_1
            ├── sample_2
            ├── sample_3
            ├── ...
            └── sample_N
```

Where `dataset` could be a sequencing RunID (e.g. `AGRF_CAGRF21077456_HFY5VDRXY`) and `sample_*` could be sample ULNs (e.g. `21_01176`).
FASTQ files for a sample are placed under it's corresponding sample directory and should use `_R1.fastq.gz` and `_R2.fastq.gz` as the filename suffix.
This is most easily done by creating symbolic links to wherever the data resides on the HPC (most likely under `/cancer/storage/SAGC/fastq/`).
Ultimately, the filename paths have the following structure:

```bash
qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}_R1.fastq.gz
qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}_R2.fastq.gz
```

If multiple replicates exist for a sample (e.g. due to multiple lanes), the workflow will aggregate the resulting BAM files into a single BAM per sample.

# Snakemake One-Time Setup

```bash
SNAKEMAKE_VERSION='6.7.0'

# Create an empty environment
conda create \
  --yes \
  --name "snakemake_v${SNAKEMAKE_VERSION}"

# Activate the new, empty environment
conda activate \
  "snakemake_v${SNAKEMAKE_VERSION}"

# Install mamba - a faster/better version of the conda executable
conda install \
  --yes \
  --channel conda-forge \
  mamba

# Install Snakemake
mamba install \
  --yes \
  --channel bioconda \
  --channel conda-forge \
  snakemake=${SNAKEMAKE_VERSION}
```
