# SGC-QC
QC pipeline for SGC GWAS data

Assumption: cohort level files are/can be named in convention cohort_PHENO_stratum for processing throughout this (i.e use Jake's cohort names). If true everything below will run automatically - otherwise needs tweaking.

1. Run qc_setup
2. Run qc_munge
3. Run METAL_example_single.metal (one script per meta-analysis)


Code not for portal: Post-processing: sgc_pipeline contains everything needed to regenerate excel doc + downstream analyses.
