# simulation data & experiment data of FSSEM
we give the analytic code of real data experiment data in paper. 
Run `run_FSSEM.R` to do FSSEM analysis of real data GSE33356

### intermediate data GE.rds and SNP.rds
we save intermediate data of FSSEM as `./inst/GE_NSCLC.rds`(Gene expression matrix) and `./inst/SNP_NSCLC.rds`(eQTL quantative matrix), which are processed GSE33356 data.


### simulation data
If you are using LSF supercomputer, you can use `run_SIM.R` to build simulation and plot simulationb result.