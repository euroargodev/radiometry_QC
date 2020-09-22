library(ggplot2)
library(ncdf4)
library(stringr)

corr_list = read.csv("~/Documents/radiometry_QC/float_list.t", sep=" ", header=T)
index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)
path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
ident = strsplit(files,"/") #separate the different roots of the files paths
ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
dac = ident[,1] #retrieve the DAC of all profiles as a vector
wod = ident[,2] #retrieve the WMO of all profiles as a vector
prof_id = ident[,4] #retrieve all profiles  name as a vector
variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file


### Get the file locations
WMO_list = corr_list$WMO[which(corr_list$corrected==1)]

file_list = rep(NA, length(WMO_list))
for (i in 1:length(WMO_list)) {
    file_A = files[wod == WMO_list[i]][1]
    file_A_sep = unlist(strsplit(file_A, "/"))
    str_sub(file_A_sep[4],2,2) = "D"
    
    file_D_sep = rep(NA, 5)
    file_D_sep[c(1,2,3,5)] = file_A_sep
    file_D_sep[4] = "RADM/RADM_profiles"
    file_D = paste0(file_D_sep, collapse = "/")
    file_list[i] = file_D
}
file_list = paste0(path_to_netcdf, file_list)

PARAMETER = "DOWN_IRRADIANCE490"

drift_C = rep(NA, length(WMO_list))
drift_Q = rep(NA, length(WMO_list))
for (i in 1:length(WMO_list)) {
    fnc = nc_open(file_list[i])
    
    parameters = ncvar_get(fnc, "STATION_PARAMETERS")
    id_param_arr = which(parameters == str_pad(PARAMETER, 64, side="right"), arr.ind=TRUE)
    if (length(id_param_arr)==2) id_prof=id_param_arr[2] else id_prof=id_param_arr[1,2]
    if (length(id_param_arr)==2) id_param=id_param_arr[1] else id_param=id_param_arr[1,1]

    #n_param = fnc$dim$N_PARAM$len
    
    sci_comment = ncvar_get(fnc, "SCIENTIFIC_CALIB_COEFFICIENT", start=c(1,id_param,1,id_prof), count=c(256,1,1,1))	
    
    nc_close(fnc)
    
    coeffs = unlist(strsplit(sci_comment, ","))
    id_C = grep("C", coeffs)
    id_Q = grep("Q", coeffs)
    
    if (length(id_C) != 0) {
        C = as.numeric(unlist(str_split(coeffs[id_C], "="))[2])
        drift_C[i] = C
    }
    if (length(id_Q) != 0) {
        Q = as.numeric(unlist(str_split(coeffs[id_Q], "="))[2])
        drift_Q[i] = Q
    }
}

ramps = drift_C[which(is.na(drift_Q) & drift_C!=0)]

p <- ggplot() +
    geom_histogram( aes(x = ramps, y=1e-7*length(ramps)*..density..), fill="#69b3a2", binwidth=1e-7 , color="#FFFFFF", boundary=T  ) +
    #geom_label( aes(x=250, y=6, label="corrected"), color="#69b3a2") +
    labs(x="dE_dark/dtime", y="Number of floats")

