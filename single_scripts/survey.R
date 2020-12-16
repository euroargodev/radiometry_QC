library(stringr)

source("~/Documents/radiometry/possol.R")

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


is_coriolis = (dac=="coriolis")

is_radiometry = rep(NA, length(prof_id))
for (i in 1:length(is_radiometry)) {
	is_radiometry[i] = any(variables[i][[1]] == "DOWNWELLING_PAR")
}

month_list = as.numeric(str_sub(prof_date,5,6))
day_list = as.numeric(str_sub(prof_date,7,8))
hour_list = as.numeric(str_sub(prof_date,9,10))
minute_list = as.numeric(str_sub(prof_date,11,12))
second_list = as.numeric(str_sub(prof_date,13,14))
tu_list = (hour_list + minute_list/60 + second_list/3600)

valid = which(!is.na(prof_date) & !is.na(lat) & !is.na(lon))

M = mapply(possol.vec, month=month_list[valid], jday=day_list[valid], tu=tu_list[valid], xlon=lon[valid], xlat=lat[valid])

solar_elev = rep(NA, length(prof_id))
solar_elev[valid] = unlist(M, use.names=F)

is_night = (solar_elev < -5)

drift_start = 20140801000000 # august 1 2014
is_drift = (prof_date > drift_start)

floats = unique(wod[which(is_coriolis & is_radiometry)])
num_prof = rep(NA, length(floats))
num_night = rep(NA, length(floats))
num_drift = rep(NA, length(floats))
is_baffin = rep(NA, length(floats))


for (i in 1:length(floats)) {
	num_prof[i] = length(which( (wod == floats[i]) & is_radiometry ))
	num_night[i] = length(which( (wod == floats[i]) & is_radiometry & is_night ))
	
	num_drift[i] = length(which( (wod == floats[i]) & is_radiometry & is_drift ))
	first_prof = which(wod == floats[i])[1]
	is_baffin[i] = (lat[first_prof] > 65) & (lon[first_prof] > -85) & (lon[first_prof] < -45)
}

num_drift[which(is_baffin)] = 0

source("~/Documents/radiometry/single_scripts/count_profiles_corrected.R")

match_WMO = match(corr_list$WMO, floats)
corr_list$num_night = num_night[match_WMO]
corr_list$num_drift = num_drift[match_WMO]
corr_list$num_prof = num_prof[match_WMO]

cat(
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof)),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof)),
"\n")
cat(
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & corr_list$corrected==1)),
"\n")
cat(
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "PEEK" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "PEEK" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "ALU" & corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
length(which(corr_list$material == "ALU" & !corr_list$num_night>0 & !corr_list$num_drift>0.8*corr_list$num_prof & !is.na(corr_list$corrected))),
"\n")
