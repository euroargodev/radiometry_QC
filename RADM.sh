#!/bin/bash

usage() { 
  echo "Usage: $0 -W <WMO_number> [-m <n_cores>] [-pz]
Do '$0 -h' for help" 1>&2
  exit 1 
}
helprint() {
  echo "
################################################################################

RADM is a tool to create and write Delayed Mode (DM) of radiometry for argo 
floats.

Usage: $0 -W <WMO_number> [-m <n_cores>] [-pz]

### Options

-W <WMO_number> : Do the delayed mode on all profiles of a float identified with
                  its 7 digits WMO number.
[-m <n_cores>] : Define the number of cores to use for parallelization. By 
                 default RADM uses all available cores.
[-p] : Plot mode, use after the DM is finished to compare before/after.
[-z] : Zoom on pressure between [0,250] decibars in plot mode.

################################################################################
" 1>&2
  exit 0
}

WMO=NA
multi_core=NA
plot_mode=FALSE
zoom_plot=FALSE

while getopts W:m:ph option
do
case "${option}"
in
W) WMO=${OPTARG};;
m) multi_core=${OPTARG};;
p) plot_mode=TRUE;;
z) zoom_plot=TRUE;;
h) helprint;;
*) usage;;
esac
done

Rscript ~/Documents/radiometry/start_RADM.R $WMO $multi_core $plot_mode $zoom_plot
