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

source("~/Documents/radiometry/possol.R")
source("~/Documents/radiometry/sensor_temp.R")
#source("~/Documents/cornec_chla_qc/chl_bbp_ttt/DM_writing/write_DM.R")

n_cores = detectCores()

#path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
path_to_netcdf = "/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)
index_greylist = read.csv("~/Documents/radiometry/ar_greylist.txt", sep=",")

files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
ident = strsplit(files,"/") #separate the different roots of the files paths
ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
dac = ident[,1] #retrieve the DAC of all profiles as a vector
wod = ident[,2] #retrieve the WMO of all profiles as a vector
prof_id = ident[,4] #retrieve all profiles  name as a vector
variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file
#variables = strsplit(variables," ") #separate the different available variables of each profile
lat = index_ifremer$latitude #retrieve the latitude of all profiles as a vector
lon = index_ifremer$longitude #retrieve the longitude of all profiles as a vector
prof_date = index_ifremer$date #retrieve the date of all profiles as a vector


#WMO = "6901525"
#WMO = "6901524"
##WMO = "6902827" # no nights ?
#WMO = "6901494" #not enough drift points ?
#WMO = "6901576"
#WMO = "6902735"
WMO = "6901473" # again drift data missing, maybe other sensor issues
#WMO = "6901474" 
#WMO = "6901495" # no drift data at all (code 290)
#WMO = "6901584"
#WMO = "6901658" # Xing works very well, new method doesn't because light penetrates very deep
#WMO = "6902547"
#WMO = "6902742"
#WMO = "6902828" # Xing works well, new method not quite on PAR
#WMO = "6902879" # no drift data (code 290)
#WMO = "6902906" # no drift data (code 290)
#WMO = "6903551" # drift bizarre, très peu de variation de Ts, bad data ?
#WMO = "7900561" # both methods work very well
#WMO = "6901492" # both good
#WMO = "6903025" # Xing great, new method fails like 6901658 because of deep light gradients



subset = which( wod==WMO & substr(prof_id,14,14)!="D" & !is.na(prof_date) & grepl("DOWNWELLING_PAR", variables))
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
day = which(solar_elev >= -5)

profile_night = profile_list[night]
files_night = files_list[night]
date_night = as.Date(as.character(prof_date_list[night]), format="%Y%m%d%H%M%S", tz="UTC")
date_night = julian(date_night, origin=as.Date("1950-01-01", tz="UTC"))

profile_day = profile_list[day]
files_day = files_list[day]
date_day = as.Date(as.character(prof_date_list[day]), format="%Y%m%d%H%M%S", tz="UTC")
date_day = julian(date_day, origin=as.Date("1950-01-01", tz="UTC"))

date_list = as.Date(as.character(prof_date_list), format="%Y%m%d%H%M%S", tz="UTC")
date_list = julian(date_list, origin=as.Date("1950-01-01", tz="UTC"))

greylist_starts = as.Date(as.character(index_greylist$START_DATE), format="%Y%m%d", tz="UTC")
index_greylist$START_JULD = julian(greylist_starts, origin=as.Date("1950-01-01", tz="UTC"))

greylist_ends = as.Date(as.character(index_greylist$END_DATE), format="%Y%m%d", tz="UTC")
index_greylist$END_JULD = julian(greylist_ends, origin=as.Date("1950-01-01", tz="UTC"))
index_greylist$END_JULD[which(is.na(index_greylist$END_JULD))] = Inf

PARAM_NAMES = c("DOWNWELLING_PAR", "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")

is_greylisted <- function(julian_day, WMO, PARAMETER_NAME) {
	relevant_lines = which(index_greylist$PLATFORM_CODE==WMO & index_greylist$PARAMETER_NAME==PARAMETER_NAME)
	greylisted_bool = any(index_greylist$START_JULD[relevant_lines] <= julian_day &
						  index_greylist$END_JULD[relevant_lines] >= julian_day)
	return(greylisted_bool)
}

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

	return(list("PARAM"=PARAM, "Ts"=fitted_Ts, "PRES"=PRES_B, "id_prof"=id_prof_B, "n_levels"=n_levels_B))
}


get_profile_match <- function(file_name, param_name, PROFILE_DATE, method="night", drift_A=0, drift_C=0, drift_Q=0) {

	match = get_Ts_match(path_to_netcdf=path_to_netcdf, file_name=file_name, PARAM_NAME=param_name)

	match_not_na = which(!is.na(match$PARAM) & !is.na(match$Ts))
	
	### correct for drift	
	param_undrifted = as.numeric( match$PARAM - drift_A - drift_C * PROFILE_DATE 
									- drift_Q * PROFILE_DATE^2)

	if (method == "night") {	

		MATCH = param_undrifted[match_not_na]
		MATCH_Ts = match$Ts[match_not_na]
		MATCH_date = rep(PROFILE_DATE, length(match_not_na))
		MATCH_name = rep(param_name, length(match_not_na))

		return(list("MATCH"=MATCH, "MATCH_Ts"=MATCH_Ts, "MATCH_date"=MATCH_date, "MATCH_name"=MATCH_name))

	} 

	if (method == "day") {
		
		if (length(match_not_na) > 4) {	

			### select dark from lilliefors test following Organelli et al
			match_param = match$PARAM[match_not_na]
			lillie_pval = rep(NA, length(match_param))
			for (j in 1:(length(match_param)-4)) {
				lillie_pval[j] = lillie.test(match_param[j:length(match_param)])$p.value
			}
			signif = (abs(lillie_pval) > 0.01)

			if (any(signif, na.rm=T)) {
				j_dark = which(signif)[1] - 1
				subsel_dark = match_not_na[j_dark:length(match_not_na)]
	
				MATCH = param_undrifted[subsel_dark]
				MATCH_Ts = match$Ts[subsel_dark]
				MATCH_date = rep(PROFILE_DATE, length(subsel_dark))
				MATCH_name = rep(param_name, length(subsel_dark))
				MATCH_darkmed = median(MATCH)
				return(list("MATCH"=MATCH, "MATCH_Ts"=MATCH_Ts, "MATCH_date"=MATCH_date, "MATCH_name"=MATCH_name,
							"MATCH_darkmed"=MATCH_darkmed))
			}
		}
		return(list("MATCH"=NULL, "MATCH_Ts"=NULL, "MATCH_date"=NULL, "MATCH_name"=NULL, "MATCH_darkmed"=NULL))
	}

	if (method == "drift") {

		match_param = match$PARAM[match_not_na]
		match_Ts = match$Ts[match_not_na]
		match_pres = match$PRES[match_not_na]
		
		valid_for_drift = which(match_pres > 240 & match_Ts < 10)
		
		if (length(valid_for_drift) == 0) {
			return(list("MATCH"=NULL, "MATCH_Ts"=NULL, "MATCH_date"=NULL, "MATCH_name"=NULL))
		}		

		MATCH = match_param[valid_for_drift]
		MATCH_Ts = match_Ts[valid_for_drift]
		MATCH_date = rep(PROFILE_DATE, length(valid_for_drift))
		MATCH_name = rep(param_name, length(valid_for_drift))

		return(list("MATCH"=MATCH, "MATCH_Ts"=MATCH_Ts, "MATCH_date"=MATCH_date, "MATCH_name"=MATCH_name))
	}	
}

####################################
### Start of drift considerations
####################################

traj_name_B = paste(path_to_netcdf, dac[subset[1]], "/", WMO, "/", WMO, "_BRtraj.nc", sep="")
traj_name_C = paste(path_to_netcdf, dac[subset[1]], "/", WMO, "/", WMO, "_Rtraj.nc", sep="")

fnc_B = nc_open(traj_name_B)
fnc_C = nc_open(traj_name_C)

temp = ncvar_get(fnc_C, "TEMP")
juld_B = ncvar_get(fnc_B, "JULD")
juld_C = ncvar_get(fnc_C, "JULD")
code = ncvar_get(fnc_B, "MEASUREMENT_CODE")

drift_C = which(!is.na(temp) & code==290) #290 may not be the only correct code


DRIFT = NULL
DRIFT_Ts = NULL
DRIFT_date = NULL
DRIFT_name = NULL
for (param_name in PARAM_NAMES) {
	irr = ncvar_get(fnc_B, param_name)

	drift_B = which(!is.na(irr) & code==290) #290 may not be the only correct code

	drift_match = rep(NA, length(drift_B))
	for (i in 1:length(drift_B)) {
		match_dist = abs(juld_B[drift_B[i]] - juld_C[drift_C])
		drift_match[i] = min(which( match_dist == min(match_dist) ))
	}

	DRIFT = c(DRIFT, irr[drift_B])	
	DRIFT_Ts = c(DRIFT_Ts, temp[drift_C[drift_match]])
	DRIFT_date = c(DRIFT_date, juld_B[drift_B])
	DRIFT_name = c(DRIFT_name, rep(param_name, length(drift_B)))	
}
nc_close(fnc_B)
nc_close(fnc_C)

drift_dataf = data.frame("PARAM"=DRIFT, "PARAM_Ts"=DRIFT_Ts, "PARAM_date"=DRIFT_date, "PARAM_name"=DRIFT_name)

#files_drift_PARALLEL = rep(files_list, 4)
#date_drift_PARALLEL = rep(date_list, 4)
#PARAM_NAMES_drift_PARALLEL = rep(PARAM_NAMES, each=length(files_list))

#DRIFT_MATCH = mcmapply(get_profile_match, file_name=files_drift_PARALLEL, param_name=PARAM_NAMES_drift_PARALLEL,
#						PROFILE_DATE=date_drift_PARALLEL, mc.cores=n_cores, USE.NAMES=FALSE, 
#						MoreArgs=list(method="drift"))
#DRIFT_dataf_all = data.frame("PARAM"=unlist(DRIFT_MATCH[1,]), 
#						"PARAM_Ts"=unlist(DRIFT_MATCH[2,]), 
#					   	"PARAM_date"=unlist(DRIFT_MATCH[3,]), 
#						"PARAM_name"=unlist(DRIFT_MATCH[4,]))


### remove drift outliers
is_drift_outlier = rep(NA, length(drift_dataf$PARAM))
#is_DRIFT_outlier = rep(NA, length(drift_dataf$PARAM))
for (param_name in PARAM_NAMES) {
	param_axis = which(drift_dataf$PARAM_name == param_name)

	quartiles = quantile(drift_dataf$PARAM[param_axis], probs=c(0.25, 0.75))
	margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
	outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)

	is_drift_outlier[param_axis] = (drift_dataf$PARAM[param_axis] < outliers_lim[1] | drift_dataf$PARAM[param_axis] > outliers_lim[2])

	### Alternative method
	#param_axis = which(DRIFT_dataf$PARAM_name == param_name)

	#quartiles = quantile(DRIFT_dataf$PARAM[param_axis], probs=c(0.25, 0.75))
	#margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
	#outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)

	#is_DRIFT_outlier[param_axis] = (DRIFT_dataf$PARAM[param_axis] < outliers_lim[1] | DRIFT_dataf$PARAM[param_axis] > outliers_lim[2])
}

drift_dataf$is_drift_outlier = is_drift_outlier
#DRIFT_dataf$is_drift_outlier = is_DRIFT_outlier

drift_dataf$PARAM_date_squared = (drift_dataf$PARAM_date)^2

drift_dataf$is_greylisted = mapply(is_greylisted, julian_day=drift_dataf$PARAM_date,
								PARAMETER_NAME=as.character(drift_dataf$PARAM_name),
								MoreArgs=list(WMO=WMO))

do_quadratic_fit = c(F, F, F, F)

A_axis_drift = rep(NA, 4)
B_axis_drift = rep(NA, 4)
C_axis_drift = rep(NA, 4)
Q_axis_drift = rep(0, 4)
drift_dataf_5C = drift_dataf
fitted_coeff_drift = NULL
#A_axis_DRIFT = rep(NA, 4)
#B_axis_DRIFT = rep(NA, 4)
#C_axis_DRIFT = rep(NA, 4)
#DRIFT_dataf_5C = DRIFT_dataf
#fitted_coeff_DRIFT = NULL
for (i in 1:4) {
	subset_PAR = which(drift_dataf$PARAM_name == PARAM_NAMES[i])
	subset_fit = which(drift_dataf$PARAM_name == PARAM_NAMES[i] & !drift_dataf$is_greylisted
						& !drift_dataf$is_drift_outlier)
	
	if (do_quadratic_fit[i]) {
		fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date + PARAM_date_squared, data=drift_dataf, subset=subset_fit)
		Q_axis_drift[i] = fit_ABC$coefficients[4]	
	} else {
		fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date, data=drift_dataf, subset=subset_fit)
	}

	fitted_coeff_drift[[PARAM_NAMES[i]]] = fit_ABC$coefficients
	A_axis_drift[i] = fit_ABC$coefficients[1]	
	B_axis_drift[i] = fit_ABC$coefficients[2]	
	C_axis_drift[i] = fit_ABC$coefficients[3]	

	drift_dataf_5C$PARAM[subset_PAR] = drift_dataf_5C$PARAM[subset_PAR] - B_axis_drift[i] * (drift_dataf_5C$PARAM_Ts[subset_PAR] - 5)

	### Alternative drift
	#subset_PAR = which(DRIFT_dataf$PARAM_name == PARAM_NAMES[i])
	#fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date, data=DRIFT_dataf, subset=subset_PAR)

	#fitted_coeff_DRIFT[[PARAM_NAMES[i]]] = fit_ABC$coefficients
	#A_axis_DRIFT[i] = fit_ABC$coefficients[1]	
	#B_axis_DRIFT[i] = fit_ABC$coefficients[2]	
	#C_axis_DRIFT[i] = fit_ABC$coefficients[3]	

	#DRIFT_dataf_5C$PARAM[subset_PAR] = DRIFT_dataf_5C$PARAM[subset_PAR] - B_axis_DRIFT[i] * (DRIFT_dataf_5C$PARAM_Ts[subset_PAR] - 5)
}

good_drift = which(!drift_dataf$is_greylisted & !drift_dataf$is_drift_outlier)
range_time = seq(min(c(drift_dataf$PARAM_date[good_drift], date_list)), max(c(drift_dataf$PARAM_date[good_drift], date_list)), length.out=100 )

data_fit_drift = data.frame(
	PARAM_name = rep(PARAM_NAMES, each=2),
	x = rep(range_time, 4),
	y = rep(A_axis_drift, each=2) + rep(B_axis_drift, each=2) * 5
		+ rep(C_axis_drift, each=2) * rep(range_time, 4) 
		+ rep(Q_axis_drift, each=2) * rep(range_time^2, 4) 
)
#data_fit_DRIFT = data.frame(
#	PARAM_name = rep(PARAM_NAMES, each=2),
#	x = rep(range_time, 4),
#	y = rep(A_axis_DRIFT, each=2) + rep(C_axis_DRIFT, each=2) * rep(range_time, 4) + rep(B_axis_DRIFT, each=2) * 5
#)

g4 = ggplot(na.omit(drift_dataf), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
	geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
	geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")
g5 = ggplot(na.omit(drift_dataf_5C), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
	geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
	geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
	#geom_point(data=function(x){x[x$is_drift_outlier & !x$is_greylisted, ]}, color="red", shape=4) +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")
#g6 = ggplot(na.omit(DRIFT_dataf), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
#	geom_point() +
#	scale_color_viridis() +
#	facet_wrap(~PARAM_name, scale="free_y")
#g7 = ggplot(na.omit(DRIFT_dataf_5C), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
#	geom_point() +
#	scale_color_viridis() +
#	facet_wrap(~PARAM_name, scale="free_y")
	
#g4_fit = g4 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
g5_fit = g5 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red") #+
#				geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")

#g6_fit = g6 + geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")
#g7_fit = g7 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red") +
#				geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")


###############################################
#### Start of radiometry/Ts matches extraction
###############################################

files_night_PARALLEL = rep(files_night, 4)
date_night_PARALLEL = rep(date_night, 4)
PARAM_NAMES_night_PARALLEL = rep(PARAM_NAMES, each=length(files_night))
A_axis_drift_night_PARALLEL = rep(A_axis_drift, each=length(files_night))
C_axis_drift_night_PARALLEL = rep(C_axis_drift, each=length(files_night))
Q_axis_drift_night_PARALLEL = rep(Q_axis_drift, each=length(files_night))

files_day_PARALLEL = rep(files_day, 4)
date_day_PARALLEL = rep(date_day, 4)
PARAM_NAMES_day_PARALLEL = rep(PARAM_NAMES, each=length(files_day))
A_axis_drift_day_PARALLEL = rep(A_axis_drift, each=length(files_day))
C_axis_drift_day_PARALLEL = rep(C_axis_drift, each=length(files_day))
Q_axis_drift_day_PARALLEL = rep(Q_axis_drift, each=length(files_day))

NIGHT_MATCH = mcmapply(get_profile_match, file_name=files_night_PARALLEL, param_name=PARAM_NAMES_night_PARALLEL,
						PROFILE_DATE=date_night_PARALLEL, drift_A=A_axis_drift_night_PARALLEL,
						drift_C=C_axis_drift_night_PARALLEL, drift_Q=Q_axis_drift_night_PARALLEL,
						mc.cores=n_cores, USE.NAMES=FALSE)
DAY_MATCH = mcmapply(get_profile_match, file_name=files_day_PARALLEL, param_name=PARAM_NAMES_day_PARALLEL,
						PROFILE_DATE=date_day_PARALLEL, drift_A=A_axis_drift_day_PARALLEL,
						drift_C=C_axis_drift_day_PARALLEL, drift_Q=Q_axis_drift_day_PARALLEL,
						mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(method="day"))

### find dark median outliers in day profiles
all_param_names = unlist(lapply(DAY_MATCH[4,], `[[`, 1)) # extract the first element of each profile
is_dark_outlier = rep(NA, length(all_param_names))
all_dark_median = unlist(DAY_MATCH[5,])
for (param_name in PARAM_NAMES) {
	param_axis = which(all_param_names == param_name)

	quartiles = quantile(all_dark_median[param_axis], probs=c(0.25, 0.75))
	margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
	outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)

	is_dark_outlier[param_axis] = (all_dark_median[param_axis] < outliers_lim[1] | all_dark_median[param_axis] > outliers_lim[2])
}
all_lengths = mapply(length, DAY_MATCH[1,])
all_lengths = all_lengths[which(all_lengths != 0)]

is_dark_outlier_full = rep(is_dark_outlier, times=all_lengths)


PAR_dataf = data.frame("PARAM"=unlist(NIGHT_MATCH[1,]), "PARAM_Ts"=unlist(NIGHT_MATCH[2,]), 
					   "PARAM_date"=unlist(NIGHT_MATCH[3,]), "PARAM_name"=unlist(NIGHT_MATCH[4,]))
DAY_dataf = data.frame("PARAM"=unlist(DAY_MATCH[1,]), 
						"PARAM_Ts"=unlist(DAY_MATCH[2,]), 
					   	"PARAM_date"=unlist(DAY_MATCH[3,]), 
						"PARAM_name"=unlist(DAY_MATCH[4,]),
						"is_dark_outlier"=is_dark_outlier_full)

PAR_dataf$is_greylisted = mapply(is_greylisted, julian_day=PAR_dataf$PARAM_date,
								PARAMETER_NAME=as.character(PAR_dataf$PARAM_name),
								MoreArgs=list(WMO=WMO))
DAY_dataf$is_greylisted = mapply(is_greylisted, julian_day=DAY_dataf$PARAM_date,
								PARAMETER_NAME=as.character(DAY_dataf$PARAM_name),
								MoreArgs=list(WMO=WMO))

g1 = ggplot(na.omit(PAR_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_date, group=PARAM_name)) +
	geom_point(data=function(x){x[!x$is_greylisted, ]}) +
	#geom_point(data=function(x){x[x$is_greylisted, ]}, color="red") +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")
g1_day = ggplot(na.omit(DAY_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_date, group=PARAM_name)) +
	geom_point(data=function(x){x[!x$is_greylisted & !x$is_dark_outlier, ]}) +
	#geom_point(data=function(x){x[x$is_greylisted & !x$is_dark_outlier, ]}, color="red") +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")

###########################
### Start of Ts fitting
###########################

fitted_coeff = NULL
A_axis = rep(NA, 4)
B_axis = rep(NA, 4)
fitted_coeff_day = NULL
A_axis_day = rep(NA, 4)
B_axis_day = rep(NA, 4)
for (i in 1:length(PARAM_NAMES)) {
	subset_PAR = which(PAR_dataf$PARAM_name==PARAM_NAMES[i] & !is.na(PAR_dataf$PARAM_Ts)
						& !PAR_dataf$is_greylisted)

	fit_AB = lm(PARAM ~ PARAM_Ts, data=PAR_dataf, subset=subset_PAR) 

	fitted_coeff[[PARAM_NAMES[i]]] = fit_AB$coefficients
	A_axis[i] = fit_AB$coefficients[[1]]
	B_axis[i] = fit_AB$coefficients[[2]]

	
	subset_DAY = which(DAY_dataf$PARAM_name==PARAM_NAMES[i] & !is.na(DAY_dataf$PARAM_Ts)
						& !DAY_dataf$is_greylisted & !DAY_dataf$is_dark_outlier)

	fit_param = DAY_dataf$PARAM[subset_DAY]
	fit_Ts = DAY_dataf$PARAM_Ts[subset_DAY]
	
	fit_param = fit_param[order(fit_Ts)]
	fit_Ts = fit_Ts[order(fit_Ts)]
	
	run_min = running(fit_param, fun=min, width=51, pad=TRUE)
	
	data_run = data.frame("run_min"=run_min, "fit_Ts"=fit_Ts)
	fit_AB_day = lm(run_min ~ fit_Ts, data=data_run)

	fitted_coeff_day[[PARAM_NAMES[i]]] = fit_AB_day$coefficients
	A_axis_day[i] = fit_AB_day$coefficients[[1]]
	B_axis_day[i] = fit_AB_day$coefficients[[2]]
	
}

Ts_range = range(c(PAR_dataf$PARAM_Ts, DAY_dataf$PARAM_Ts), na.rm=T)

data_fit = data.frame(
	PARAM_name = rep(PARAM_NAMES, each=2),
	x = rep(Ts_range, 4),
	y = rep(A_axis, each=2) + rep(B_axis, each=2) * rep(Ts_range, 4)
)
data_fit_day = data.frame(
	PARAM_name = rep(PARAM_NAMES, each=2),
	x = rep(Ts_range, 4),
	y = rep(A_axis_day, each=2) + rep(B_axis_day, each=2) * rep(Ts_range, 4)
)

g2 = g1 + geom_line(data=data_fit, mapping=aes(x=x,y=y), color="red")
g3 = g2 + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="black")

g2_day = g1_day + geom_line(data=data_fit, mapping=aes(x=x,y=y), color="red")
g3_day = g2_day + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="black")

save_plots=TRUE
if (save_plots) {
	plot_name=paste0(WMO, "_regr_day.png")
	png(filename=plot_name, width=600, height=600)
	plot(g3_day)
	dev.off()

	plot_name=paste0(WMO, "_regr_night.png")
	png(filename=plot_name, width=600, height=600)
	plot(g3)
	dev.off()

	plot_name=paste0(WMO, "_regr_drift.png")
	png(filename=plot_name, width=600, height=600)
	plot(g5_fit)
	dev.off()

	plot_name=paste0(WMO, "_regr_drift_uncorrected.png")
	png(filename=plot_name, width=600, height=600)
	plot(g4)
	dev.off()
}

#all_match_380 = mcmapply(get_Ts_match, file_name=files_list, mc.cores=n_cores, SIMPLIFY=FALSE,
#							MoreArgs=list(path_to_netcdf=path_to_netcdf, PARAM_NAME="DOWN_IRRADIANCE380"))

#plot_corr_380 <- function(n) {
#
#	corr = all_match_380[[n]]$PARAM - ( A_axis[2] + B_axis[2] * all_match_380[[n]]$Ts )
#
#	dataf3 = data.frame(IRR380=all_match_380[[n]]$PARAM, PRES=all_match_380[[n]]$PRES, corr=corr, Ts=all_match_380[[n]]$Ts)
#
#	g3 = ggplot(na.omit(dataf3), aes(x=IRR380, y=PRES)) +
#		geom_point() +
#		scale_x_continuous(trans="log10") +
#		scale_y_continuous(trans="reverse") +
#		geom_path(mapping=aes(x=corr, y=PRES), color="red")
#	
#	g3_1 = ggplot(na.omit(dataf3), aes(x=IRR380, y=PRES)) +
#		geom_point() +
#		scale_x_continuous(limits=c(-1e-4, 1e-4)) +
#		scale_y_continuous(trans="reverse") +
#		geom_path(mapping=aes(x=corr, y=PRES), color="red")
#
#	g4 = ggplot(na.omit(dataf3), aes(x=Ts, y=PRES)) +
#		geom_point() + 
#		scale_y_continuous(trans="reverse")
#
#	grid.arrange(g3_1, g3, g4, nrow=1)
#
#}

####################################
### Apply correction to netcdf
####################################

date_update = Sys.time()
DATE = stri_datetime_format(date_update, format="uuuuMMddHHmmss", tz="UTC")

corr_file <- function(file_name, PROFILE_DATE, path_to_netcdf, use_day=FALSE) {

	path_sep = unlist(strsplit(file_name, "/"))
	if (use_day) {
		file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], "radiometry_xing_day", path_sep[4], sep="/")
	} else {
		file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], "radiometry_xing", path_sep[4], sep="/")
	}
	full_file_name_out = paste(path_to_netcdf, file_name_out, sep="")
	full_file_name_in = paste(path_to_netcdf, file_name, sep="")
	
	system2("cp", c(full_file_name_in, full_file_name_out))

	for (i in 1:length(PARAM_NAMES)) {
		
		matchup = get_Ts_match(file_name=file_name, path_to_netcdf=path_to_netcdf, PARAM_NAME=PARAM_NAMES[i])
		
		param_undrifted = as.numeric( matchup$PARAM - A_axis_drift[i] - C_axis_drift[i] * PROFILE_DATE
										- Q_axis_drift[i] * PROFILE_DATE^2 )
		
		if (use_day) {
			corr = param_undrifted - ( A_axis_day[i] + B_axis_day[i] * matchup$Ts )
		} else {
			corr = param_undrifted - ( A_axis[i] + B_axis[i] * matchup$Ts )
		}
		
		corr_error = 0.2*corr #TODO

		corr_qc = rep("1", length(corr)) #TODO
		corr_qc[which(is.na(corr))] = "4"
		corr_qc = paste(corr_qc, collapse="")				

		fnc = nc_open(full_file_name_out, write=TRUE)
		
		ncvar_put(fnc, paste(PARAM_NAMES[i],"_ADJUSTED",sep=""), corr, start=c(1, matchup$id_prof), count=c(matchup$n_levels, 1))
		
		nc_close(fnc)
	}
	return(0)	
		
}

#corr_file(files_list[1], path_to_netcdf)

#corr_all = mcmapply(corr_file, file_name=files_list, PROFILE_DATE=date_list, mc.cores=n_cores, SIMPLIFY=FALSE,
#							MoreArgs=list(path_to_netcdf=path_to_netcdf, use_day=FALSE))



