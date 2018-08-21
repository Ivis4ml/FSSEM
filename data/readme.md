# Simulation data & experiment data of FSSEM
we give the analytic code of real data experiment data in paper. 
Run `run_FSSEM.R` to do FSSEM analysis of real data GSE33356

### Intermediate data for NSCLC real data analysis
we save intermediate data of FSSEM as `./inst/GE_NSCLC.rds`(Gene expression matrix) and `./inst/SNP_NSCLC.rds`(eQTL quantative matrix), which are processed GSE33356 data.


### Simulation data
If you are using LSF server, you can use `./inst/run_SIM.R` to build simulation jobs to submit to LSF server and then plot simulation result.