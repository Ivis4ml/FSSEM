source("grandnet.R")            ### generate simulated DAG/DCG for FSSEM
source("datastd.R")             ### data centralization 
source("ridge.R")               ### ridge regression for FSSEM's initial
source("fssem.R")               ### FSSEM algorithm & cross-validation related functions
source("util.R")                ### utility funcitons
source("stabsel.R")             ### stability selection

### compressed version solver of FSSEM integrated in solver.min.R. 
### For convinience and mobility, directly source solver.min.R to run FSSEM.
### source("solver.min.R") 