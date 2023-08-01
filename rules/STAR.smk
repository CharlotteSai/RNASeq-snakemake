import math

R1_FILENAME_SUFFIX = "_R1.fastq.gz"
R2_FILENAME_SUFFIX = "_R2.fastq.gz"
SE_FILENAME_SUFFIX = '.fastq.gz'

STAR_INDEXING_MEM_GB = 48

localrules:
	STAR_GenomeLength,
	STAR_NumberOfReferences,

STAR_INDEX_FILES = [
  'chrLength.txt',
  'chrNameLength.txt',
  'chrName.txt',
  'chrStart.txt',
  'SA',
  'SAindex',
  'Genome',
  'genomeParameters.txt',
]

rule STAR_indexing:
	input:
		ref                = "{fasta}",
		gtf                = "{fasta}.gtf",
		NumberOfReferences = "{fasta}.gz_STAR_index/NumberOfReferences",
		GenomeLength       = "{fasta}.gz_STAR_index/GenomeLength"
	output:
		expand("{{fasta}}.gz_STAR_index/{file}", file=STAR_INDEX_FILES)
	threads:
		24,
	resources:
		threads = 24,
		mem_gb  = STAR_INDEXING_MEM_GB,
		time_hr = 24,
	params:
		tmp_dir                = os.environ["LOCAL_SCRATCH"],
		limitGenomeGenerateRAM = STAR_INDEXING_MEM_GB * 1000000000
	conda:
		"../envs/STAR.yaml",
	benchmark:
		repeat("benchmarks/STAR_indexing/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		TMP_OUT="$(mktemp --directory --dry-run {params.tmp_dir}/tmp.XXXXXXXXXX_STARtmp)"
		TMP_GENOME_OUT="$(mktemp --directory {params.tmp_dir}/tmp.XXXXXXXXXX_STAR_index)"

		function clean_up {{
		  # Perform program exit housekeeping
		  if [ ! -z "${{TMP_OUT}}" ] && [ -e "${{TMP_OUT}}" ]; then
		    rm -rf "${{TMP_OUT}}" "${{TMP_GENOME_OUT}}"
		  fi
		  trap 0  # reset to default action
		  exit
		}}
		trap clean_up 0 1 2 3 15 #see 'man signal' for descriptions http://man7.org/linux/man-pages/man7/signal.7.html
 
		genomeChrBinNbits=$(echo "v=l($(cat {input.GenomeLength})/$(cat {input.NumberOfReferences})) / l(2); scale=0; def=18; if (def<v) def else v/1" | bc -l)
		genomeSAindexNbases=$(echo "v=(l($(cat {input.GenomeLength}))/2-1) / l(2); def=14; if(def<v) def else v/1; scale=0;" | bc -l)

		STAR \
		  --runMode genomeGenerate \
		  --outTmpDir ${{TMP_OUT}} \
		  --genomeDir ${{TMP_GENOME_OUT}} \
		  --genomeChrBinNbits ${{genomeChrBinNbits}} \
		  --genomeSAindexNbases ${{genomeSAindexNbases}} \
		  --genomeFastaFiles {input.ref} \
		  --sjdbGTFfile {input.gtf} \
		  --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} \
		  --runThreadN {threads}

		mv ${{TMP_GENOME_OUT}}/* {wildcards.fasta}.gz_STAR_index/
		"""

# TODO: potentially mark output as temp() so than once these files are merged, we don't duplicate storage requirements
# TODO: using genomeLoad = LoadAndKeep generates this error:
#   Shared memory error: 4, errno: Invalid argument(22)
#   EXITING because of FATAL ERROR: problems with shared memory: error from shmget() or shm_open().
#   SOLUTION: check shared memory settings as explained in STAR manual, OR run STAR with --genomeLoad NoSharedMemory to avoid using shared memory

checkpoint STAR_mapping_PE:
	input:
		R1    = "qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R1_FILENAME_SUFFIX,
		R2    = "qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R2_FILENAME_SUFFIX,
		index = expand("references/{{fasta}}_STAR_index/{file}", file=STAR_INDEX_FILES),
	output:
		bam              = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.Aligned.sortedByCoord.out.bam",
		splice_junctions = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.SJ.out.tab",
		flag             = touch("mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.done"),
	log:
		main     = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.Log.out",
		progress = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.Log.progress.out",
		final    = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.Log.final.out",
		std      = "mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/{sample}/{replicate}.Log.std.out",
	wildcard_constraints:
		sample  = "[^\/]+",
		dataset = "[^\/]+",
	resources:
		threads = 20,
		mem_gb  = lambda wildcards, attempt, input: math.ceil( (sum(os.path.getsize(f) for f in input['index'] if os.path.isfile(f)) / 1024**3) + (0.30 * 20) * attempt ),
		time_hr = lambda wildcards, attempt, input: max(math.ceil(sum(os.path.getsize(f) for f in [input['R1'],input['R2']] if os.path.isfile(f)) / 1024**3 / 30 / 60), 1) * attempt,
	threads:
		20,
	params:
		tmp_dir                     = os.environ["LOCAL_SCRATCH"],
		outBAMsortingThreadN        = lambda wildcards, threads: max(1, threads - 6),
		limitIObufferSize_I         = 50 * 1000000,
		limitIObufferSize_O         = 100 * 1000000,
		limitBAMsortRAM             = lambda wildcards, resources: resources['mem_gb'] * 1024**3 * 0.9,
		outFileNamePrefix           = lambda wildcards, output: output['bam'].replace('Aligned.sortedByCoord.out.bam', ''),
		outFilterMultimapNmax       = 5,
		alignEndsType               = 'Local',
		outFilterScoreMinOverLread  = 0.5,
		outFilterMatchNminOverLread = 0.5,
		outFilterMismatchNoverLmax  = 0.3,
		alignIntronMax              = 10000,
		alignMatesGapMax            = 10000,
		genomeLoad                  = 'NoSharedMemory',
	conda:
		"../envs/STAR.yaml",
	benchmark:
		repeat("benchmarks/STAR_mapping_PE/{fasta}{dataset}{sample}{replicate}.txt", N_BENCHMARKS),
	shell:
		"""
		TMP_OUT="$(mktemp --directory --dry-run {params.tmp_dir}/tmp.XXXXXXXXXX_STARtmp)"

		function clean_up {{
		  # Perform program exit housekeeping
		  if [ ! -z "${{TMP_OUT}}" ] && [ -e "${{TMP_OUT}}" ]; then
		    rm -rf "${{TMP_OUT}}"
		  fi
		  trap 0  # reset to default action
		  exit
		}}
		trap clean_up 0 1 2 3 15 #see 'man signal' for descriptions http://man7.org/linux/man-pages/man7/signal.7.html

		STAR \
		  --runMode alignReads \
		  --outStd BAM_SortedByCoordinate \
		  --outBAMcompression 0 \
		  --outTmpDir ${{TMP_OUT}}/ \
		  --genomeDir references/{wildcards.fasta}_STAR_index \
		  --genomeLoad {params.genomeLoad} \
		  --runThreadN {threads} \
		  --outBAMsortingThreadN {params.outBAMsortingThreadN} \
		  --limitIObufferSize {params.limitIObufferSize_I} {params.limitIObufferSize_O} \
		  --limitBAMsortRAM {params.limitBAMsortRAM} \
		  --readFilesIn {input.R1} {input.R2} \
		  --readFilesCommand pigz -dcp2 \
		  --outFileNamePrefix {params.outFileNamePrefix} \
		  --outSAMtype BAM SortedByCoordinate \
		  --outSAMstrandField intronMotif \
		  --outSAMattributes All \
		  --outSAMattrRGline ID:{output[0]} PL:Illumina PU:Unknown LB:Unknown SM:{wildcards.sample} \
		  --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} --outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax} \
		  --outSJfilterOverhangMin 35 20 20 20 \
		  --outSJfilterCountTotalMin 10 3 3 3 \
		  --outSJfilterCountUniqueMin 5 1 1 1 \
		  --alignIntronMax {params.alignIntronMax} \
		  --alignMatesGapMax {params.alignMatesGapMax} \
		| samtools calmd -b --threads {threads} - references/{wildcards.fasta} \
		> {output.bam}

		samtools quickcheck -v {output.bam}
		"""

#def multiqc_input(wildcards):
#	INPUTS = glob.glob('mapped_reads/' + wildcards.fasta + '/star/transcriptomic/PE/' + wildcards.dataset + '/*/*.Log.final.out')
#
#	return(INPUTS)


rule multiqc_STAR_PE:
	input:
		# BUG: this runs at the start of DAG generation and returns empty as mapping hasn't been concluded yet
		lambda wildcards: glob.glob('mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}/*/*.Log.final.out'.format(fasta=wildcards.fasta, dataset=wildcards.dataset))
#		expand(
#			'mapped_reads/{{fasta}}/star/transcriptomic/PE/{{dataset}}/{sample}/{replicate}.Log.final.out',
#		),
#		'mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}',
	output:
		html = 'mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}_multiqc_report.html',
		data = directory('mapped_reads/{fasta}/star/transcriptomic/PE/{dataset}_multiqc_report_data'),
	params:
		max_table_rows = 10000,
	conda:
		'../envs/multiqc.yaml',
	benchmark:
		repeat('benchmarks/multiqc_STAR_PE/{fasta}/{dataset}/benchmark.txt', N_BENCHMARKS),
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 ** attempt,
		threads = 1,
	threads:
		1,
	shell:
		"""
		multiqc {input} \
		  --module star \
		  --force \
		  --filename {output.html} \
		  --cl_config '{{max_table_rows: {params.max_table_rows}}}'
		"""

rule STAR_mapping_SE:
	input:  
		SE    = "qc_reads/transcriptomic/SE/{dataset}/{sample}/{replicate}" + SE_FILENAME_SUFFIX,
		index = expand("references/{{fasta}}_STAR_index/{file}", file=STAR_INDEX_FILES),
	output: 
		bam	             = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.Aligned.sortedByCoord.out.bam",
		splice_junctions = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.SJ.out.tab"
	log:
		main     = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.Log.out",
		progress = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.Log.progress.out",
		final    = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.Log.final.out",
		std      = "mapped_reads/{fasta}/star/transcriptomic/SE/{dataset}/{sample}/{replicate}.Log.std.out",
	wildcard_constraints:
		sample  = "[^\/]+",
		dataset = "[^\/]+",
	resources:
		threads = 20,
		mem_gb  = lambda wildcards, input: math.ceil( (sum(os.path.getsize(f) for f in input['index'] if os.path.isfile(f)) / 1024**3) + (0.15 * MAX_THREADS) ) + 5,
		time_hr = lambda wildcards, attempt, input: max(math.ceil(sum(os.path.getsize(f) for f in [input['SE']] if os.path.isfile(f)) / 1024**3 / 30 / 60), 15) * 10 * attempt
	threads:
		20,
	params:
		tmp_dir                     = os.environ["LOCAL_SCRATCH"],
		outBAMsortingThreadN        = 14,
		limitIObufferSize           = 150 * 1000000,
		limitBAMsortRAM             = 190 * 1000000000,
		outFileNamePrefix           = lambda wildcards, output: output['bam'].replace('Aligned.sortedByCoord.out.bam', ''),
		outFilterMultimapNmax       = 5,
		alignEndsType               = 'Local',
		outFilterMismatchNoverLmax  = 0.3,
		outFilterMatchNminOverLread = 1.00 - 0.3,
		alignIntronMax              = 10000,
		alignMatesGapMax            = 10000,
		genomeLoad                  = 'NoSharedMemory'
	benchmark:
		repeat("benchmarks/STAR_mapping_SE/{fasta}{dataset}{sample}{replicate}.txt", N_BENCHMARKS),
	shell:
		"""
		TMP_OUT="$(mktemp --directory --dry-run {params.tmp_dir}/tmp.XXXXXXXXXX_STARtmp)"

		function clean_up {{
		  # Perform program exit housekeeping
		  if [ ! -z "${{TMP_OUT}}" ] && [ -e "${{TMP_OUT}}" ]; then
		    rm -rf "${{TMP_OUT}}"
		  fi
		  trap 0  # reset to default action
		  exit
		}}
		trap clean_up 0 1 2 3 15 #see 'man signal' for descriptions http://man7.org/linux/man-pages/man7/signal.7.html

		STAR \
		  --runMode alignReads \
		  --outStd BAM_SortedByCoordinate \
		  --outBAMcompression 0 \
		  --outTmpDir ${{TMP_OUT}}/ \
		  --genomeDir references/{wildcards.fasta}_STAR_index \
		  --genomeLoad {params.genomeLoad} \
		  --runThreadN {threads} \
		  --outBAMsortingThreadN {params.outBAMsortingThreadN} \
		  --limitIObufferSize {params.limitIObufferSize} \
		  --limitBAMsortRAM {params.limitBAMsortRAM} \
		  --readFilesIn {input.SE} \
		  --readFilesCommand pigz -dcp2 \
		  --outFileNamePrefix {params.outFileNamePrefix} \
		  --outSAMtype BAM SortedByCoordinate \
		  --outFilterMultimapScoreRange 0 \
		  --outFilterMultimapNmax {params.outFilterMultimapNmax} \
		  --outFilterMismatchNoverLmax {params.outFilterMismatchNoverLmax} \
		  --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
		  --outSJfilterOverhangMin 35 20 20 20 \
		  --outSJfilterCountTotalMin 10 3 3 3 \
		  --outSJfilterCountUniqueMin 5 1 1 1 \
		  --alignEndsType {params.alignEndsType} \
		  --alignSoftClipAtReferenceEnds No \
		  --outSAMstrandField intronMotif \
		  --outSAMattributes All \
		  --alignIntronMax {params.alignIntronMax} \
		  --alignMatesGapMax {params.alignMatesGapMax} \
		  --outSAMattrRGline ID:{output[0]} PL:Illumina PU:Unknown LB:Unknown SM:{wildcards.sample} \
		| samtools calmd -b --threads {threads} - references/{wildcards.fasta} \
		> {output.bam}

		samtools quickcheck -v {output.bam}
		"""

rule STAR_NumberOfReferences:
	input:
		"{fasta}.fai",
	output:
		"{fasta}.gz_STAR_index/NumberOfReferences",
	resources:
		threads = 1,
		mem_gb  = 1,
		time_hr = 1,
	shell:
		"""
		wc -l < {input} > {output}
		"""

rule STAR_GenomeLength:
	input:
		"{fasta}.fai",
	output:
		"{fasta}.gz_STAR_index/GenomeLength",
	resources:
		threads = 1,
		mem_gb  = 1,
		time_hr = 1,
	conda:
		"../envs/mawk.yaml",
	shell:
		"""
		mawk '{{tot+=$2}}END{{printf "%d", tot}}' {input} > {output}
		"""
