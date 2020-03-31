library(ncdf4)
library(stringr)
library(parallel)

source("~/Documents/radiometry/possol.R")

#path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
path_to_netcdf = "/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)

files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
ident = strsplit(files,"/") #separate the different roots of the files paths
ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
dac = ident[,1] #retrieve the DAC of all profiles as a vector
wod = ident[,2] #retrieve the WMO of all profiles as a vector
prof_id = ident[,4] #retrieve all profiles  name as a vector
variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file
variables = strsplit(variables," ") #separate the different available variables of each profile
lat = index_ifremer$latitude #retrieve the latitude of all profiles as a vector
lon = index_ifremer$longitude #retrieve the longitude of all profiles as a vector
prof_date = index_ifremer$date #retrieve the date of all profiles as a vector


WMO = "6901524"
#WMO = "6902827"

profile_list = substr(prof_id[which(wod==WMO)], 3, 14)
files_list = files[which(wod==WMO)]
lat_list = lat[which(wod==WMO)]
lon_list = lon[which(wod==WMO)]
prof_date_list = prof_date[which(wod==WMO)]

month_list = as.numeric(str_sub(prof_date_list,5,6))
day_list = as.numeric(str_sub(prof_date_list,7,8))
hour_list = as.numeric(str_sub(prof_date_list,9,10))
minute_list = as.numeric(str_sub(prof_date_list,11,12))
second_list = as.numeric(str_sub(prof_date_list,13,14))
tu_list = (hour_list + minute_list/60 + second_list/3600)

M = mapply(possol.vec, month=month_list, jday=day_list, tu=tu_list, xlon=lon_list, xlat=lat_list)

solar_elev = unlist(M, use.names=F)

night = which(solar_elev < -5)

profile_night = profile_list[night]
files_night = files_list[night]

get_xing_AB <- function(path_to_netcdf, file_name, PARAM) {
	
	file_B = paste(path_to_netcdf, file_name, sep="")
	
	# get core file name
	path_split = unlist(strsplit(file_name, "/"))
	path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
	filenc_name_C = paste("?",substring(path_split[4], 3),sep="")
	file_C = paste(path_to_netcdf, path_to_profile, "/", filenc_name_C, sep="") 
	file_C = system2("ls", file_C, stdout=TRUE) # identify R or D file 
	if (length(file_C)==2) { # if both R and D files exist
		file_C = file_C[1] # use the D file which is first in alphabetical order
	}
		
	filenc_B = nc_open(file_B)
	filenc_C = nc_open(file_C)
	

	nc_close(filenc_B)
	nc_close(filenc_C)
}

ab = get_xing_AB(path_to_netcdf, files_night[1], "DOWNWELLING_PAR") 
