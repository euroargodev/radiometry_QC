############################################################################################################
### Functions to extract matches between Ts (sensor temperature) and radiometry parameters in 
### argo profiles
############################################################################################################


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
    n_prof_B = filenc_B$dim$N_PROF$len
    
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
    
    return(list("PARAM"=PARAM, "Ts"=fitted_Ts, "PRES"=PRES_B, "id_prof"=id_prof_B, "n_prof"=n_prof_B, "n_levels"=n_levels_B,
                "PARAM_QC"=PARAM_QC))
}

get_profile_match <- function(file_name, param_name, path_to_netcdf, PROFILE_DATE, method="night", drift_A=0, drift_C=0, drift_Q=0) {
    
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
                            "MATCH_darkmed"=MATCH_darkmed, "MATCH_which"=subsel_dark))
            }
        }
        return(list("MATCH"=NULL, "MATCH_Ts"=NULL, "MATCH_date"=NULL, "MATCH_name"=NULL, "MATCH_darkmed"=NULL,
                    "MATCH_which"=NULL))
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
    
    if (method == "drift_xing") {
        
        match_param = match$PARAM[match_not_na]
        match_Ts = match$Ts[match_not_na]
        match_pres = match$PRES[match_not_na]
        
        run_sd = as.vector(running(match_param, fun=sd, width=10, pad=T))
        
        #lim_param = list("DOWN_IRRADIANCE380" = 5e-5,
        #                "DOWN_IRRADIANCE412" = 4e-5,
        #                 "DOWN_IRRADIANCE490" = 2e-5,
        #                 "DOWNWELLING_PAR" = 5e-2) # test 1
        #lim_param = list("DOWN_IRRADIANCE380" = 1e-4,
        #                 "DOWN_IRRADIANCE412" = 1e-4,
        #                 "DOWN_IRRADIANCE490" = 5e-5,
        #                 "DOWNWELLING_PAR" = 5e-2) # test 2
        lim_param = list("DOWN_IRRADIANCE380" = 5e-5,
                         "DOWN_IRRADIANCE412" = 5e-5,
                         "DOWN_IRRADIANCE490" = 5e-5,
                         "DOWNWELLING_PAR" = 5e-2) # just orders of magnitude management
        #lim_param = list("DOWN_IRRADIANCE380" = 3e-5,
        #                 "DOWN_IRRADIANCE412" = 7e-5,
        #                 "DOWN_IRRADIANCE490" = 7e-5,
        #                 "DOWNWELLING_PAR" = 5e-2) # accounting for typical proportionality (6901658)
        #lim_param = list("DOWN_IRRADIANCE380" = 1e-4,
        #                 "DOWN_IRRADIANCE412" = 1e-3,
        #                 "DOWN_IRRADIANCE490" = 1e-2,
        #                 "DOWNWELLING_PAR" = 5e-2) # accounting for typical proportionality (7900561)
        
        lim = max(which(run_sd >= lim_param[[param_name]]))
        
        if (lim == length(match_param) | lim == -Inf) {
            return(list("MATCH"=NULL, "MATCH_Ts"=NULL, "MATCH_date"=NULL, "MATCH_name"=NULL))
        }	
        
        MATCH = match_param[(lim+1):length(match_param)]
        MATCH_Ts = match_Ts[(lim+1):length(match_param)]
        MATCH_date = rep(PROFILE_DATE, length(match_param)-lim)
        MATCH_name = rep(param_name, length(match_param)-lim)
        
        return(list("MATCH"=MATCH, "MATCH_Ts"=MATCH_Ts, "MATCH_date"=MATCH_date, "MATCH_name"=MATCH_name))
        
    }
    
    if (method=="test") {
        MATCH = match$PARAM
        MATCH_Ts = match$Ts
        MATCH_date = rep(PROFILE_DATE, length(match$PARAM))
        MATCH_name = rep(param_name, length(match$PARAM))
        
        return(list("MATCH"=MATCH, "MATCH_Ts"=MATCH_Ts, "MATCH_date"=MATCH_date, "MATCH_name"=MATCH_name))
        
    }
    
}
