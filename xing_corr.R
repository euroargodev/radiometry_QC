library(ncdf4)
library(stringr)
library(stringi)
library(parallel)
library(ggplot2)
library(viridis)
library(gridExtra)
library(gtools)
library(nortest)

source("~/Documents/radiometry/possol.R")
source("~/Documents/radiometry/sensor_temp.R")
#source("~/Documents/cornec_chla_qc/chl_bbp_ttt/DM_writing/write_DM.R")

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
day = which(solar_elev >= -5)

profile_night = profile_list[night]
files_night = files_list[night]
date_night = as.Date(as.character(prof_date_list[night]), format="%Y%m%d%H%M%S", origin="1950-01-01", tz="UTC")

profile_day = profile_list[day]
files_day = files_list[day]
date_day = as.Date(as.character(prof_date_list[day]), format="%Y%m%d%H%M%S", origin="1950-01-01", tz="UTC")


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

PARAM_NAMES = c("DOWNWELLING_PAR", "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")

all_PARAM_match = NULL

PAR = NULL
PAR_Ts = NULL
PAR_date = NULL
PAR_name = NULL
DAY = NULL
DAY_Ts = NULL
DAY_date = NULL
DAY_name = NULL

for (param_name in PARAM_NAMES) {
	for (i in 1:length(files_night)) {
		match = get_Ts_match(path_to_netcdf, files_night[i], param_name)

		PAR = c(PAR, match$PARAM)
		PAR_Ts = c(PAR_Ts, match$Ts)

		PAR_date_new = rep(date_night[i], length(match$PARAM))
		PAR_date_new[which(is.na(match$PARAM))] = NA
		PAR_date = c(PAR_date, PAR_date_new)
		
		PAR_name_new = rep(param_name, length(match$PARAM))
		PAR_name_new[which(is.na(match$PARAM))] = NA
		PAR_name = c(PAR_name, PAR_name_new)
	}
	#PAR_date = as.Date(PAR_date, origin="1970-01-01")
	for (i in 1:length(files_day)) {
		match = get_Ts_match(path_to_netcdf, files_day[i], param_name)
		
		match_not_na = which(!is.na(match$PARAM))
		if (length(match_not_na) > 4) {
			
			### select dark from lillifors test following Organelli et al
			match_param = match$PARAM[match_not_na]
			lillie_pval = rep(NA, length(match_param))
			for (j in 1:(length(match_param)-4)) {
				lillie_pval[j] = lillie.test(match_param[j:length(match_param)])$p.value
			}
			signif = (abs(lillie_pval) > 0.01)
			if (any(signif, na.rm=T)) {
				j_dark = which(signif)[1] - 1
				subsel_dark = match_not_na[j_dark:length(match_not_na)]

				DAY = c(DAY, match$PARAM[subsel_dark])
				DAY_Ts = c(DAY_Ts, match$Ts[subsel_dark])
	
				DAY_date_new = rep(date_day[i], length(subsel_dark))
				#DAY_date_new[which(is.na(match$PARAM))] = NA
				DAY_date = c(DAY_date, DAY_date_new)
			
				DAY_name_new = rep(param_name, length(subsel_dark))
				#DAY_name_new[which(is.na(match$PARAM))] = NA
				DAY_name = c(DAY_name, DAY_name_new)
			}
		}
	}
}
PAR_dataf = data.frame("PARAM"=PAR, "PARAM_Ts"=PAR_Ts, "PARAM_date"=PAR_date, "PARAM_name"=PAR_name)
DAY_dataf = data.frame("PARAM"=DAY, "PARAM_Ts"=DAY_Ts, "PARAM_date"=DAY_date, "PARAM_name"=DAY_name)

g1 = ggplot(na.omit(PAR_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_date, group=PARAM_name)) +
	geom_point() +
	scale_color_viridis() +
	#scale_y_continuous(limits=c(-1e-4, -2e-5)) +
	facet_wrap(~PARAM_name, scale="free_y")
g1_day = ggplot(na.omit(DAY_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_date, group=PARAM_name)) +
	geom_point() +
	scale_color_viridis() +
	#scale_y_continuous(limits=c(-1e-4, -2e-5)) +
	facet_wrap(~PARAM_name, scale="free_y")

fitted_coeff = NULL
A_axis = rep(NA, 4)
B_axis = rep(NA, 4)
fitted_coeff_day = NULL
A_axis_day = rep(NA, 4)
B_axis_day = rep(NA, 4)
for (i in 1:length(PARAM_NAMES)) {
	subset_PAR = which(PAR_dataf$PARAM_name==PARAM_NAMES[i] & !is.na(PAR_dataf$PARAM_Ts))

	fit_AB = lm(PARAM ~ PARAM_Ts, data=PAR_dataf, subset=subset_PAR) 

	fitted_coeff[[PARAM_NAMES[i]]] = fit_AB$coefficients
	A_axis[i] = fit_AB$coefficients[[1]]
	B_axis[i] = fit_AB$coefficients[[2]]

	
	fit_param = DAY_dataf$PARAM[which(DAY_dataf$PARAM_name == PARAM_NAMES[i])]
	fit_Ts = DAY_dataf$PARAM_Ts[which(DAY_dataf$PARAM_name == PARAM_NAMES[i])]
	
	fit_param = fit_param[order(fit_Ts)]
	fit_Ts = fit_Ts[order(fit_Ts)]
	
	run_min = running(fit_param, fun=min, width=51, pad=TRUE)
	
	data_run = data.frame("run_min"=run_min, "fit_Ts"=fit_Ts)
	fit_AB_day = lm(run_min ~ fit_Ts, data=data_run)

	fitted_coeff_day[[PARAM_NAMES[i]]] = fit_AB_day$coefficients
	A_axis_day[i] = fit_AB_day$coefficients[[1]]
	B_axis_day[i] = fit_AB_day$coefficients[[2]]
	
}

Ts_range = range(PAR_dataf$PARAM_Ts, na.rm=T)

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
g2_1 = g2 + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="green")

g2_day = g1_day + geom_line(data=data_fit, mapping=aes(x=x,y=y), color="red")
g2_1_day = g2_day + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="green")


n_cores = detectCores()
all_match_380 = mcmapply(get_Ts_match, file_name=files_list, mc.cores=n_cores, SIMPLIFY=FALSE,
							MoreArgs=list(path_to_netcdf=path_to_netcdf, PARAM_NAME="DOWN_IRRADIANCE380"))

plot_corr_380 <- function(n) {

	corr = all_match_380[[n]]$PARAM - ( A_axis[2] + B_axis[2] * all_match_380[[n]]$Ts )

	dataf3 = data.frame(IRR380=all_match_380[[n]]$PARAM, PRES=all_match_380[[n]]$PRES, corr=corr, Ts=all_match_380[[n]]$Ts)

	g3 = ggplot(na.omit(dataf3), aes(x=IRR380, y=PRES)) +
		geom_point() +
		scale_x_continuous(trans="log10") +
		scale_y_continuous(trans="reverse") +
		geom_path(mapping=aes(x=corr, y=PRES), color="red")
	
	g3_1 = ggplot(na.omit(dataf3), aes(x=IRR380, y=PRES)) +
		geom_point() +
		scale_x_continuous(limits=c(-1e-4, 1e-4)) +
		scale_y_continuous(trans="reverse") +
		geom_path(mapping=aes(x=corr, y=PRES), color="red")

	g4 = ggplot(na.omit(dataf3), aes(x=Ts, y=PRES)) +
		geom_point() + 
		scale_y_continuous(trans="reverse")

	grid.arrange(g3_1, g3, g4, nrow=1)

}

####################################
### Apply correction to netcdf
####################################

date_update = Sys.time()
DATE = stri_datetime_format(date_update, format="uuuuMMddHHmmss", tz="UTC")

corr_file <- function(file_name, path_to_netcdf, use_day=FALSE) {

	path_sep = unlist(strsplit(file_name, "/"))
	if (use_day) {
		file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], "radiometry_xing_day", path_sep[4], sep="/")
	} else {
		file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], "radiometry_xing", path_sep[4], sep="/")
	}
	full_file_name_out = paste(path_to_netcdf, file_name_out, sep="")
	full_file_name_in = paste(path_to_netcdf, file_name, sep="")
	
	system2("cp", c(full_file_name_in, full_file_name_out))

	for (param_name in PARAM_NAMES) {
		
		matchup = get_Ts_match(file_name=file_name, path_to_netcdf=path_to_netcdf, PARAM_NAME=param_name)
		
		if (use_day) {
			corr = matchup$PARAM - ( fitted_coeff_day[[param_name]][1] + fitted_coeff_day[[param_name]][2] * matchup$Ts )
		} else {
			corr = matchup$PARAM - ( fitted_coeff[[param_name]][1] + fitted_coeff[[param_name]][2] * matchup$Ts )
		}
		
		corr_error = 0.2*corr #TODO

		corr_qc = rep("1", length(corr)) #TODO
		corr_qc[which(is.na(corr))] = "4"
		corr_qc = paste(corr_qc, collapse="")				

		fnc = nc_open(full_file_name_out, write=TRUE)
		
		ncvar_put(fnc, paste(param_name,"_ADJUSTED",sep=""), corr, start=c(1, matchup$id_prof), count=c(matchup$n_levels, 1))
		
		nc_close(fnc)
	}
	return(0)	
		
}

#corr_file(files_list[1], path_to_netcdf)

#corr_all = mcmapply(corr_file, file_name=files_list, mc.cores=n_cores, SIMPLIFY=FALSE,
#							MoreArgs=list(path_to_netcdf=path_to_netcdf, use_day=TRUE))

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


PAR = NULL
PAR_Ts = NULL
PAR_date = NULL
PAR_name = NULL
for (param_name in PARAM_NAMES) {
	irr = ncvar_get(fnc_B, param_name)

	drift_B = which(!is.na(irr) & code==290) #290 may not be the only correct code

	drift_match = rep(NA, length(drift_B))
	for (i in 1:length(drift_B)) {
		match_dist = abs(juld_B[drift_B[i]] - juld_C[drift_C])
		drift_match[i] = min(which( match_dist == min(match_dist) ))
	}

	#k = 0.19/60 # s^-1
	#delta_t = 54 # s
	#lag = (delta_t + 1/k) / (3600*24) # d # lag for Ts when trailing a linear Tw

	#Ts = rep(NA, length(drift_C))
	#Ts[1] = temp[drift_C[1]]
	#for (i in 2:length(drift_C)) {
	#	Ts[i] = temp[drift_C[i]] - lag * (temp[drift_C[i]] - temp[drift_C[i-1]]) / (juld_C[drift_C[i]] - juld_C[drift_C[i-1]])
	#} # we assume a long time between drift measurement and linear rate of change between them
	
	PAR = c(PAR, irr[drift_B])	
	PAR_Ts = c(PAR_Ts, temp[drift_C[drift_match]])
	PAR_date = c(PAR_date, juld_B[drift_B])
	PAR_name = c(PAR_name, rep(param_name, length(drift_B)))	
}
nc_close(fnc_B)
nc_close(fnc_C)

drift_dataf = data.frame("PARAM"=PAR, "PARAM_Ts"=PAR_Ts, "PARAM_date"=PAR_date, "PARAM_name"=PAR_name)

A_axis_drift = rep(NA, 4)
B_axis_drift = rep(NA, 4)
C_axis_drift = rep(NA, 4)
drift_dataf_5C = drift_dataf
for (i in 1:4) {
	subset_PAR = which(drift_dataf$PARAM_name == PARAM_NAMES[i])
	fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date, data=drift_dataf, subset=subset_PAR)
	A_axis_drift[i] = fit_ABC$coefficients[1]	
	B_axis_drift[i] = fit_ABC$coefficients[2]	
	C_axis_drift[i] = fit_ABC$coefficients[3]	

	drift_dataf_5C$PARAM[subset_PAR] = drift_dataf_5C$PARAM[subset_PAR] - B_axis_drift[i] * (drift_dataf_5C$PARAM_Ts[subset_PAR] - 5)
}

range_time = range(drift_dataf$PARAM_date)
data_fit_drift = data.frame(
	PARAM_name = rep(PARAM_NAMES, each=2),
	x = rep(range_time, 4),
	y = rep(A_axis_drift, each=2) + rep(C_axis_drift, each=2) * rep(range_time, 4) + rep(B_axis_drift, each=2) * 5
)

g4 = ggplot(na.omit(drift_dataf), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
	geom_point() +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")
g5 = ggplot(na.omit(drift_dataf_5C), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
	geom_point() +
	scale_color_viridis() +
	facet_wrap(~PARAM_name, scale="free_y")
	
g4_fit = g4 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
g5_fit = g5 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")


####################################
#### End of drift considerations
####################################


