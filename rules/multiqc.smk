rule multiqc:
	input:
		compose_input_for_multiqc_raw_r1
	output:
		html = INPUT_PROJECT_DIR + '/QC-results/multiqc_report.html',
		data = directory(INPUT_PROJECT_DIR + '/QC-results/multiqc_report_data'),
	params:
		max_table_rows = 10000,
	conda:
		"envs/qc.yml",
	benchmark:
		repeat(BENCHMARKS_STORAGE + PIPELINE_NAME + '/' + PIPELINE_VERSION + "/multiqc_raw_r1" + INPUT_PROJECT_DIR + "/benchmark.txt", N_BENCHMARKS),
	threads:
		1,
	shell:
		"""
		multiqc {input} --force --filename {output.html} --cl_config '{{max_table_rows: {params.max_table_rows}}}'
		"""
