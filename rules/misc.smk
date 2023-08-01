rule bai_index_bam:
	input:
		"{bam}",
	output:
		"{bam}.bai",
	wildcard_constraints:
		bam = "^.+?\.bam",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 4 * attempt,
		threads = 12,
	threads:
		12,
	conda:
		"../envs/samtools.yaml",
	benchmark:
		repeat("benchmarks/bai_index_bam/{bam}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools index -b -@ {threads} {input}
		"""

rule csi_index_bam:
	input:
		"{bam}",
	output:
		"{bam}.csi",
	wildcard_constraints:
		bam = "^.+?\.bam",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 4 * attempt,
		threads = 12,
	conda:
		"../envs/samtools.yaml",
	threads:
		12,
	benchmark:
		repeat("benchmarks/csi_index_bam/{bam}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools index -c -@ {threads} {input}
		"""

rule uncompress_fasta:
	input:
		"{fasta}.gz",
	output:
		"{fasta}",
	wildcard_constraints:
		fasta = ".+?(fasta|fas|fa)",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 4,
	threads:
		4,
	conda:
		"../envs/pigz.yaml",
	benchmark:
		repeat("benchmarks/uncompress_fasta/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		pigz -dcp{threads} {input} > {output}
		"""

rule create_fasta_fai:
	input:
		"{fasta}",
	output:
		"{fasta}.fai",
	wildcard_constraints:
		fasta = ".+?(fasta|fas|fa)",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 2,
	threads:
		2,
	conda:
		"../envs/samtools.yaml",
	benchmark:
		repeat("benchmarks/create_fasta_fai/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {input}
		"""

rule create_fasta_gzi:
	input:
		"{fasta}",
	output:
		"{fasta}.gzi",
	wildcard_constraints:
		fasta = ".+?(fasta|fas|fa)\.gz",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 2,
	threads:
		2,
	benchmark:
		repeat("benchmarks/create_fasta_gzi/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools faidx {input}
		"""

rule make_chrom_size:
	input:
		"{fasta}.fai",
	output:
		"{fasta}.chrom.sizes",
	wildcard_constraints:
		fasta = ".+?(fasta|fas|fa)",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 1,
	benchmark:
		repeat("benchmarks/make_chrom_size/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		cut -f1,2 {input} > {output}
		"""

rule fasta_gaps_as_GFF3:
	input:
		ref = "{fasta}",
	output:
		gff3 = "{fasta}.gaps.gff3.gz",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 4 * attempt,
		threads = 12,
	threads:
		12,
	conda:
		"../envs/fasta_gaps_as_GFF3.yml",
	benchmark:
		repeat("benchmarks/fasta_gaps_as_GFF3/{fasta}.txt", N_BENCHMARKS),
	shell:
		"""
		pigz -dcp2 {input.ref} \
		  | fasta_formatter \
		  | parallel --will-cite --pipe --max-lines 2 --keep-order --max-procs {threads} ./scripts/fasta2gff3Gaps.pl \
		  | sed '1i ##gff-version 3\n' \
		  | bgzip --threads {threads} \
		  > {output.gff3}
		"""

rule tabix_index_gff3_tbi:
	input:
		gff3 = "{prefix}.gff3.gz",
	output:
		index = "{prefix}.gff3.gz.tbi",
	conda:
		"../envs/tabix_index_gff3_tbi.yml",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 1,
	benchmark:
		repeat("benchmarks/tabix_index_gff3_tbi/{prefix}.txt", N_BENCHMARKS),
	shell:
		"""
		tabix --preset gff {input.gff3}
		"""

rule tabix_index_gff3_csi:
	input:
		gff3 = "{prefix}.gff3.gz",
	output:
		index = "{prefix}.gff3.gz.csi",
	conda:
		"../envs/tabix_index_gff3_tbi.yml",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 1,
	benchmark:
		repeat("benchmarks/tabix_index_gff3_csi/{prefix}.txt", N_BENCHMARKS),
	shell:
		"""
		tabix --csi --preset gff {input.gff3}
		"""

rule bam_to_cram:
	input:
		bam = "by_seq/{prefix}/{ref}/{file}.bam",
		ref = "references/{ref}",
	output:
		"by_seq/{prefix}/{ref}/{file}.cram",
	params:
		seqs_per_slice = 1000,
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 8 * attempt,
		threads = 12,
	threads:
		12,
	benchmark:
		repeat("benchmarks/bam_to_cram/{prefix}{ref}{file}.txt", N_BENCHMARKS),
	shell:
		"""
		scramble -t {threads} \
		  -r "{input.ref}" \
		  -I bam \
		  -O cram \
		  -s {params.seqs_per_slice} \
		  "{input.bam}" \
		  "{output}"
		"""

rule index_cram:
	input:
		"{file}.cram",
	output:
		"{file}.cram.crai",
	resources:
		time_hr = lambda wildcards, attempt: 1 ** attempt,
		mem_gb  = lambda wildcards, attempt: 1 * attempt,
		threads = 12,
	threads:
		12,
	benchmark:
		repeat("benchmarks/index_cram/{file}.txt", N_BENCHMARKS),
	shell:
		"""
		samtools index -@ {threads} {input}
		"""
