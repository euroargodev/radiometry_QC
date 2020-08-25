######################################################################
# Script to be called from bash by RADM.sh
# Starts the Delayed Mode with the adequate options
######################################################################

uf = commandArgs(trailingOnly = TRUE)

### get all arguments
WMO = uf[1]
multi_core = uf[2]

### Check conflicting options
if ( WMO=="NA" ) {
    print("Please give a WMO number as an argument (-W), or read the help (-h)")
    stop()
}


### import pathways
source("~/Documents/radiometry/pathways.R")

### Source the DM functions and subfunctions, and libraries
library(ncdf4)
library(stringr)
library(stringi)
library(parallel)
library(ggplot2)
library(plyr)
library(viridis)
library(gridExtra)
library(gtools)
library(nortest)

source(paste(path_to_source, "possol.R", sep=""))
source(paste(path_to_source, "sensor_temp.R", sep=""))
source(paste(path_to_source, "xing_corr.R", sep=""))
source(paste(path_to_source, "get_matches.R", sep=""))
source(paste(path_to_source, "write_DM.R", sep=""))
source(paste(path_to_source, "increment_N_CALIB.R", sep=""))

### import tables
index_ifremer = read.table(path_to_index_ifremer, sep=",", header = T)
index_greylist = read.csv(path_to_index_greylist, sep = ",")


### Treat optional arguments

if (multi_core == "NA") {
    n_cores = detectCores()
} else {
    n_cores = as.numeric(multi_cores)
}



### Compute and write delayed modes
exit = main_RADM(WMO=WMO, index_ifremer=index_ifremer, index_greylist=index_greylist, path_to_netcdf=path_to_netcdf, n_cores=n_cores)


