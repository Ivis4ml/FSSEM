# Simulation data & experiment data of FSSEM
we give the analytic code of real data experiment data in paper. 
Run `run_FSSEM.R` to do FSSEM analysis of real data GSE33356

## Intermediate data for NSCLC real data analysis
we save intermediate data of FSSEM as `./inst/GE_NSCLC.rds`(Gene expression matrix) and `./inst/SNP_NSCLC.rds`(eQTL quantative matrix), which are processed GSE33356 data.

### Preequirments
+ `mkdir exp`
+ `cd exp`
+ Downloading `GSE33356_family.soft.gz` from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33356
+ Unzip `inst/HG-U133_Plus_2.na36.annot.csv`
+ Downloadiing `inst/GenomeWideSNP_6.na29.annot.zip` from http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
+ Installing package `hgu133plus2hsentrezgcdf` downloaded from http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133plus2hsentrezgprobe_22.0.0.tar.gz



## Simulation data
If you are using LSF server, you can use `./inst/run_SIM.R` to build simulation jobs to submit to LSF server and then plot simulation result.