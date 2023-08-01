import os, glob

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Default filename suffixes are NOT parameterised but instead specified in:
#   compose_input_for_merge_replicates_for_sample()
#R1_FILENAME_SUFFIX = "_R1.qc.fastq.gz"
#R2_FILENAME_SUFFIX = "_R2.qc.fastq.gz"

def compose_input_for_merge_replicates_for_sample(wildcards):
	#eprint(wildcards)

	QC_SAMPLE_SUFFIX = '_R[12].fastq.gz'
	BAM_SUFFIX = '.bam'

	if wildcards.path.endswith('/SE'):
		# SE input
		QC_SAMPLE_SUFFIX = '.fq.gz'

	if wildcards.mapping_params == "star":
		# STAR transcriptomics
		BAM_SUFFIX = '.Aligned.sortedByCoord.out.bam'

	# Get a list of files from the qc_reads dir matching the location defined in wildcards
	# combined with the suffix expected for the QC files
	GLOB_PATTERN = 'qc_reads/' + wildcards.path + '/' + wildcards.dataset + '/' + wildcards.sample + '/*' + QC_SAMPLE_SUFFIX
	QC_FASTQ = glob.glob(GLOB_PATTERN)

	#eprint(GLOB_PATTERN)
	#eprint(QC_FASTQ)

	# Construct the list of BAM files by manipulating the paths of the qc_read FASTQ files
	REPLACEMENT = "mapped_reads/%(ref)s/%(mapping_params)s/\g<1>" % wildcards + BAM_SUFFIX
	regex = re.compile(r"^qc_reads/(.+?)" + QC_SAMPLE_SUFFIX.replace('.', '\.') + "$")

	INPUTS = list(set([regex.sub(REPLACEMENT, fastq) for fastq in QC_FASTQ]))
	#eprint(INPUTS)

	return(INPUTS)

rule merge_replicates_for_sample:
	input:
		compose_input_for_merge_replicates_for_sample
	output:
		bam   = "mapped_reads_merged/{ref}/{mapping_params}/{path}/{dataset}/{sample}.bam",
	wildcard_constraints:
		ref            = "[^\/]+",
		mapping_params = "[^\/]+",
		dataset        = "[^\/]+",
		sample         = "[^\/]+"
	resources:
		time_hr       = lambda wildcards, input, attempt: 1 * attempt if len(input) == 1 else 1 * attempt,
		gres_tmpfs_gb = lambda wildcards, input, attempt: 1 if len(input) == 1 else math.ceil( sum(os.path.getsize(f) for f in input if os.path.isfile(f))  / 1024**3 / 1.5 * attempt),
		mem_gb        = lambda wildcards, input, attempt: 1 * attempt if len(input) == 1 else 1 * attempt,
		threads       = 24,
	threads:
		24,
	params:
		tmp_dir           = os.environ["LOCAL_SCRATCH"],
		compression_level = 2
	conda:
		"../envs/samtools.yaml"
	priority:
		30,
	benchmark:
		repeat("benchmarks/merge_replicates_for_sample/{ref}/{mapping_params}/{path}/{dataset}/{sample}.txt", N_BENCHMARKS),
	shell:
		"""
		TMP_OUT="$(mktemp {params.tmp_dir}/tmp.XXXXXXXXXX.bam)"

		function clean_up {{
		  # Perform program exit housekeeping
		  if [ ! -z "${{TMP_OUT}}" ] && [ -e "${{TMP_OUT}}" ]; then
		    rm -f "${{TMP_OUT}}"
		  fi
		  trap 0  # reset to default action
		  exit
		}}
		trap clean_up 0 1 2 3 15 #see 'man signal' for descriptions http://man7.org/linux/man-pages/man7/signal.7.html

		INPUTS=( {input} )
		if [ "${{#INPUTS[@]}}" -eq 1 ]; then
		  # Only 1 input, so no need to merge, just symlink
		  # NOTE: Assumes the target of the symlink is not marked as temp
		  #ln --symbolic $(readlink -f {input}) {output.bam}
		  cp -a {input} {output.bam}
		else
		  samtools merge -f --threads {threads} -l {params.compression_level} ${{TMP_OUT}} {input}
		  mv "${{TMP_OUT}}" {output.bam}
		  chmod 0644 {output.bam}
		fi

		samtools quickcheck -v {output.bam}
		"""
