######################################################################
# Script to be called from bash by RADM.sh
# Starts the Delayed Mode with the adequate options
######################################################################

uf = commandArgs(trailingOnly = TRUE)

### get all arguments
WMO = uf[1]
multi_core = uf[2]
plot_mode = as.logical(uf[3])
zoom_plot = as.logical(uf[4])

### Check conflicting options
if ( WMO=="NA" ) {
    print("Please give a WMO number as an argument (-W), or read the help (-h)")
    stop()
}

cat("Importing libraries and source code...")

### import pathways
source("~/Documents/radiometry/pathways.R")

### Source the DM functions and subfunctions, and libraries
suppressPackageStartupMessages({
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
    library(tidyr)
    library(dplyr)
})
    
source(paste0(path_to_source, "possol.R"))
source(paste0(path_to_source, "sensor_temp.R"))
source(paste0(path_to_source, "main_RADM.R"))
source(paste0(path_to_source, "get_matches.R"))
source(paste0(path_to_source, "write_DM.R"))
source(paste0(path_to_source, "increment_N_CALIB.R"))
source(paste0(path_to_source, "RT_QC_radiometry_function_oao_2.R"))
source(paste0(path_to_source, "plots.R"))


cat("DONE\nImporting bio index and greylist...")

### import tables
index_ifremer = read.table(path_to_index_ifremer, sep=",", header = T)
index_greylist = read.csv(path_to_index_greylist, sep = ",")

cat("DONE\n")

### Treat optional arguments

if (multi_core == "NA") {
    n_cores = detectCores()
} else {
    n_cores = as.numeric(multi_cores)
}


if (plot_mode) {
    exit = plot_corr_wrapper(WMO, index_ifremer, path_to_netcdf, n_cores=n_cores, pres_zoom=zoom_plot)
} else {
    ### Compute and write delayed modes
    exit = main_RADM(WMO=WMO, index_ifremer=index_ifremer, index_greylist=index_greylist, path_to_netcdf=path_to_netcdf, n_cores=n_cores)
}

