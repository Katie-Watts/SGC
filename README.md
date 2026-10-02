# SGC-QC
QC pipeline for SGC GWAS data

Assumption: cohort level files are/can be named in convention cohort_PHENO_stratum for processing throughout this (i.e use Jake's cohort names).

1. Run qc_setup
2. Run qc_munge
3. Run METAL (METAL_example_single.metal - shows the parameters we want across all runs, just need to change the cohorts going into each i.e the PROCESS lines).


Code not for portal: Post-processing: sgc_pipeline contains everything needed to regenerate excel doc + downstream analyses.
