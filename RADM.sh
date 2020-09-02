#!/bin/bash

usage() { 
	echo "Usage: $0 -W <WMO_number> [-m <n_cores>]
Do '$0 -h' for help" 1>&2
	exit 1 
}
helprint() {
	echo "
#########################################################################################

RADM is a tool to create and write Delayed Mode (DM) of radiometry for argo floats.

Usage: $0 -W <WMO_number> [-m <n_cores>]

### Options

-W <WMO_number> : Do the delayed mode on all profiles of a float identified with its
                  7 digits WMO number.
[-m <n_cores>] : Define the number of cores to use for parallelization. By default RADM
                 uses all available cores.

#########################################################################################
" 1>&2
	exit 0
}

WMO=NA
multi_core=NA

while getopts W:m:h option
do
case "${option}"
in
W) WMO=${OPTARG};;
m) multi_core=${OPTARG};;
h) helprint;;
*) usage;;
esac
done

Rscript ~/Documents/radiometry/start_RADM.R $WMO $multi_core
