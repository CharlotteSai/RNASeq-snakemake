import pandas
import os

envvars:
	"LOCAL_SCRATCH",

###########################
# Project Specific Config #
###########################
configfile: "config.yaml"

REFERENCE   = config['REFERENCE']
RNA_DATASET = config['RNA_DATASET']
RNA_MAPPER  = config['RNA_MAPPER']

CHROM_IDS    = pandas.read_csv('references/' + REFERENCE + '.fai', sep="\t", names=['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH']) .NAME.tolist()
SAMPLES      = next(os.walk('./qc_reads/transcriptomic/PE/' + RNA_DATASET + '/'))[1]
N_BENCHMARKS = 1

singularity:
	"docker://continuumio/miniconda3:4.10.3p0"

include: "rules/misc.smk"
include: "rules/STAR.smk"
include: "rules/merge_samples.smk"
include: "rules/bbtools.smk"

localrules:
	all,

rule all:
	input:
		# Merged BAM files
		expand(
			'mapped_reads_merged/{ref}/{rna_mapper}/transcriptomic/PE/{dataset}/{sample}{ext}',
			ref        = REFERENCE,
			rna_mapper = RNA_MAPPER,
			dataset    = RNA_DATASET,
			sample     = SAMPLES,
			ext        = [
				'.bam',
				'.bam.csi',
			],
		),
		# MultiQC report
		expand(
			'mapped_reads/{ref}/{rna_mapper}/transcriptomic/PE/{dataset}_multiqc_report.html',
			ref        = REFERENCE,
			rna_mapper = RNA_MAPPER,
			dataset    = RNA_DATASET,
		),
