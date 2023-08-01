rule bbduk_ribo_rm:
	input:
		r1 = "qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R1_FILENAME_SUFFIX,
		r2 = "qc_reads/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R2_FILENAME_SUFFIX,
		Ref = "/data/bioinformatics/bcbio_genomes/others/rRNA_contamination/rRNA-db-contam.fasta.gz"
	output:
		r1 = "bbduk/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R1_FILENAME_SUFFIX,
		r2 = "bbduk/transcriptomic/PE/{dataset}/{sample}/{replicate}" + R2_FILENAME_SUFFIX,
	#conda:
	#	'../envs/bbtools.yaml',
	resources:
		mem_gb  = 4,
		time_hr = 1,
		threads = 2,
	params:
		bbduk_kmer = 31,
	threads:
		2,
	shell:
		"""
		/homes/charlotte.sai/bbmap/bbduk.sh -Xmx3g \
		  in={input.r1} \
		  in2={input.r2} \
		  ref={input.Ref} \
		  out={output.r1} \
		  out2={output.r2} \
		  k={params.bbduk_kmer}
		"""
