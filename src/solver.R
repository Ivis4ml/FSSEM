source("grandnet.R")            ### generate simulated DAG/DCG for FSSEM
source("datastd.R")             ### data centralization 
source("ridge.R")               ### ridge regression for FSSEM's initial
source("fssem.R")               ### FSSEM algorithm & cross-validation related functions
source("util.R")                ### utility funcitons
source("stabsel.R")             ### stability selection

########## compressed solver of FSSEM integrated in solver_min. For convinience, directly source solver_min.R
########## source("solver_min.R") 