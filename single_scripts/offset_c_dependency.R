library(ggplot2)
library(ncdf4)
library(stringr)
library(parallel)

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
WMO = "6901654"


file_A = files[wod == WMO]
file_A = file_A[which(str_sub(file_A, -4, -4) != "D")]
file_list_10 = rep(NA, length(file_A))
file_list_08 = rep(NA, length(file_A))
file_list_12 = rep(NA, length(file_A))
for (i in 1:length(file_A)) {
    
    file_A_sep = unlist(strsplit(file_A[i], "/"))
    str_sub(file_A_sep[4],2,2) = "D"
    
    file_D_sep = rep(NA, 5)
    file_D_sep[c(1,2,3,5)] = file_A_sep
    
    file_D_sep[4] = "RADM/RADM_profiles"
    file_D = paste0(file_D_sep, collapse = "/")
    file_list_10[i] = file_D
    
    file_D_sep[4] = "RADM/RADM_profiles_08cm"
    file_D = paste0(file_D_sep, collapse = "/")
    file_list_08[i] = file_D
    
    file_D_sep[4] = "RADM/RADM_profiles_12cm"
    file_D = paste0(file_D_sep, collapse = "/")
    file_list_12[i] = file_D
}
file_list_10 = paste0(path_to_netcdf, file_list_10)
file_list_08 = paste0(path_to_netcdf, file_list_08)
file_list_12 = paste0(path_to_netcdf, file_list_12)

PARAM_NAME = "DOWN_IRRADIANCE490"

get_param <- function(filename, PARAM_NAME) {
    fnc = nc_open(filename)
    
    parameters = ncvar_get(fnc, "STATION_PARAMETERS")
    id_param_arr = which(parameters == str_pad(PARAM_NAME, 64, side="right"), arr.ind=TRUE)
    if (length(id_param_arr)==2) id_prof=id_param_arr[2] else id_prof=id_param_arr[1,2]
    if (length(id_param_arr)==2) id_param=id_param_arr[1] else id_param=id_param_arr[1,1]
    
    n_level = fnc$dim$N_LEVEL$len
    
    PARAM_NAME_ADJUSTED = paste0(PARAM_NAME, "_ADJUSTED")
    
    param_adjusted = ncvar_get(fnc, PARAM_NAME_ADJUSTED, start=c(1,id_prof), count=c(n_level,1))	
    
    nc_close(fnc)
    
    return(param_adjusted)
}

n_cores = detectCores()

all_param_10 = unlist(mcmapply(get_param, filename = file_list_10, MoreArgs = list("PARAM_NAME"=PARAM_NAME), mc.cores = n_cores, USE.NAMES = F))
all_param_08 = unlist(mcmapply(get_param, filename = file_list_08, MoreArgs = list("PARAM_NAME"=PARAM_NAME), mc.cores = n_cores, USE.NAMES = F))
all_param_12 = unlist(mcmapply(get_param, filename = file_list_12, MoreArgs = list("PARAM_NAME"=PARAM_NAME), mc.cores = n_cores, USE.NAMES = F))

h_08 <- ggplot() +
    geom_histogram( aes(x = all_param_08-all_param_10, ..density..), fill="#69b3a2", color="#FFFFFF", boundary=T  )
h_12 <- ggplot() +
    geom_histogram( aes(x = all_param_12-all_param_10, ..density..), fill="#69b3a2", color="#FFFFFF", boundary=T  )

p_08 = ggplot() +
    geom_point(aes(x=all_param_10, y=all_param_08-all_param_10)) +
    scale_y_continuous(limits = c(-3e-5, -0.5e-5)) +
    scale_x_continuous(limits = c(0, 2.5e-5))
    #geom_line(aes(x=abs(all_param_10), y=abs(all_param_10)), color="red") +
    #scale_x_log10() +
    #scale_y_log10()

out_08 = which(all_param_08-all_param_10 < -1e-5)
out_08_rel = (all_param_08[out_08]-all_param_10[out_08])/all_param_10[out_08]

