library(ncdf4)
library(stringr)
library(parallel)

source("~/Documents/radiometry/possol.R")
source("~/Documents/radiometry/sensor_temp.R")

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


WMO = "6901525"
#WMO = "6901524"
#WMO = "6902827"

subset = which( wod==WMO & substr(prof_id,14,14)!="D" )
#subset = which( wod==WMO )

profile_list = substr(prof_id[subset], 3, 14)
files_list = files[subset]
lat_list = lat[subset]
lon_list = lon[subset]
prof_date_list = prof_date[subset]

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

get_Ts_match <- function(path_to_netcdf, file_name, PARAM_NAME) {
	
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
	
    parameters_B = ncvar_get(filenc_B,"STATION_PARAMETERS") 
    parameters_C = ncvar_get(filenc_C,"STATION_PARAMETERS") 
    id_param_arr_B = which(parameters_B==str_pad(PARAM_NAME, 64, side="right"), arr.ind=TRUE)
    id_param_arr_C = which(parameters_C==str_pad("TEMP", 16, side="right"), arr.ind=TRUE)
	if (length(id_param_arr_B)==2) id_prof_B=id_param_arr_B[2] else id_prof_B=id_param_arr_B[1,2]
	if (length(id_param_arr_C)==2) id_prof_C=id_param_arr_C[2] else id_prof_C=id_param_arr_C[1,2]
	n_levels_B = filenc_B$dim$N_LEVELS$len
	n_levels_C = filenc_C$dim$N_LEVELS$len # should be the same as n_levels_B
	
	PARAM = ncvar_get(filenc_B, PARAM_NAME, start=c(1,id_prof_B), count=c(n_levels_B,1))	
	TEMP = ncvar_get(filenc_C, "TEMP", start=c(1,id_prof_C), count=c(n_levels_C,1))	
	PRES_B = ncvar_get(filenc_B, "PRES", start=c(1,id_prof_B), count=c(n_levels_B,1))	
	PRES_C = ncvar_get(filenc_C, "PRES", start=c(1,id_prof_C), count=c(n_levels_C,1))	

	PARAM_QC = unlist(str_split(ncvar_get(filenc_B, paste(PARAM_NAME,"_QC",sep=""), start=c(1,id_prof_B), count=c(n_levels_B,1)), ""))
	TEMP_QC = unlist(str_split(ncvar_get(filenc_C, "TEMP_QC", start=c(1,id_prof_C), count=c(n_levels_C,1)), ""))
	PRES_B_QC = unlist(str_split(ncvar_get(filenc_C, "PRES_QC", start=c(1,id_prof_B), count=c(n_levels_C,1)), ""))
	PRES_C_QC = unlist(str_split(ncvar_get(filenc_C, "PRES_QC", start=c(1,id_prof_C), count=c(n_levels_C,1)), ""))
	
	bad_data_B = which(PARAM_QC=="3" | PARAM_QC=="4" | PRES_B_QC=="3" | PRES_B_QC=="4")
	bad_data_C = which(TEMP_QC=="3" | TEMP_QC=="4" | PRES_C_QC=="3" | PRES_C_QC=="4")
	PARAM[bad_data_B] = NA
	TEMP[bad_data_C] = NA
	PRES_B[bad_data_B] = NA
	PRES_C[bad_data_C] = NA
	
	nc_close(filenc_B)
	nc_close(filenc_C)

	fitted_Ts = sensor_temp(TEMP, PRES_C, PRES_B)

	#return(list("a"=TEMP_QC, "b"=PARAM_QC, "c"=PRES_B_QC, "d"=PRES_C_QC))
	return(list("PARAM"=PARAM, "Ts"=fitted_Ts))
}

#ab = get_Ts_match(path_to_netcdf, files_night[1], "DOWNWELLING_PAR")
#print(ab)

PAR = NULL
PAR_Ts = NULL

for (file_name in files_night) {
	match = get_Ts_match(path_to_netcdf, file_name, "DOWNWELLING_PAR")
	PAR = c(PAR, match$PARAM)
	PAR_Ts = c(PAR_Ts, match$Ts)
}

plot(PAR_Ts, PAR)
