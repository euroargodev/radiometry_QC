library(ncdf4)
library(stringr)
library(nortest)
library(parallel)

source("~/Documents/radiometry/get_matches.R")
source("~/Documents/radiometry/sensor_temp.R")

PARAM_NAMES = c("DOWNWELLING_PAR", "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")


#WMO_list = system2("ls", c("/DATA/correct_RADM/RADM_???????.zip"), stdout = T)
WMO_list = system2("ls", c("/mnt/c/DATA/correct_RADM/RADM_???????.zip"), stdout = T)
WMO_list = str_sub(WMO_list, -11, -5)
#file_list = system2("ls", c("/DATA/correct_RADM/RADM_???????/RADM_profiles/B*"), stdout = T)

corr_file <- function(file_name, WMO) {
    
    path_sep = unlist(strsplit(file_name, "/"))
    
    #file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], path_sep[4], "RADM_profiles_2", path_sep[6], sep="/")
    file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], path_sep[4], path_sep[5], path_sep[6], "RADM_profiles_2", path_sep[8], sep="/")
    full_file_name_out = paste(file_name_out, sep="")
    #full_file_name_in = file_name
    
    system2("cp", c(file_name, full_file_name_out))
 
    fnc = nc_open(full_file_name_out, readunlim=FALSE, write=TRUE)
    
    parameters = ncvar_get(fnc, "STATION_PARAMETERS")
    id_param_arr = which(parameters == str_pad("DOWNWELLING_PAR", 64, side="right"), arr.ind=TRUE)

	if (length(id_param_arr)==0) return(0)
    if (length(id_param_arr)==2) id_prof=id_param_arr[2] else id_prof=id_param_arr[1,2]
    if (length(id_param_arr)==2) id_param=id_param_arr[1] else id_param=id_param_arr[1,1]
    
    PAR_datamode = ncvar_get(fnc, "PARAMETER_DATA_MODE", start=c(id_param, id_prof), count=c(1,1))
    
    if (PAR_datamode != "D") { #PAR not in delayed mode means we did not correct this file
        nc_close(fnc)
        return(0) 
    }
    
	file_C = paste0("/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/coriolis/", WMO, "/profiles/", "?", substr(path_sep[8], 3, 16))	
    file_C = system2("ls", file_C, stdout=TRUE) # identify R or D file 
    if (length(file_C)==2) { # if both R and D files exist
       	file_C = file_C[1] # use the D file which is first in alphabetical order
    }

    for (i in 1:length(PARAM_NAMES)) {
        
        PARAM_ADJUSTED_NAME = paste0(PARAM_NAMES[i], "_ADJUSTED")
        PARAM_ADJUSTED_QC_NAME = paste0(PARAM_ADJUSTED_NAME, "_QC")
        PARAM_ADJUSTED_ERROR_NAME = paste0(PARAM_ADJUSTED_NAME, "_ERROR")
       	
		#material = "PEEK"
		material = "Aluminium"
		 
        matchup = get_Ts_match(file_name=file_name, path_to_netcdf="", PARAM_NAME=PARAM_NAMES[i], material=material, core_file_name=file_C)
        
        match_day = get_profile_match(file_name=file_name, param_name=PARAM_NAMES[i], path_to_netcdf="", 
                                      PROFILE_DATE=0, method="day", material=material, core_file_name=file_C)
        
        ### Flags
        
        corr_qc = matchup$PARAM_QC
        
        corr_qc[match_day$MATCH_which] = rep(as.character(2), length(match_day$MATCH_which)) #if which==NULL --> no change
        corr_qc[which(matchup$PARAM_QC=="3" | matchup$PARAM_QC=="4")] = "4"
        corr_qc[which(matchup$PRES_B_QC=="3" | matchup$PRES_B_QC=="4")] = "4"
        corr_qc[which(is.na(matchup$Ts) & !is.na(matchup$PARAM))] = "4"
        corr_qc_word = paste(corr_qc, collapse="")	
        
        blank_qc = paste(rep(" ", matchup$n_levels), collapse="")
        corr_qc_array = rep(blank_qc, matchup$n_prof)
        corr_qc_array[matchup$id_prof] = corr_qc_word
        
        ### Adjusted parameter
        
        #param_undrifted = as.numeric( matchup$PARAM - A_axis_drift_corr[i] - C_axis_drift_corr[i] * PROFILE_DATE
        #                              - Q_axis_drift_corr[i] * PROFILE_DATE^2 )
        
        #corr = param_undrifted - ( A_axis_corr[i] + B_axis_corr[i] * matchup$Ts )
        #corr[which(corr_qc == "4")] = NA #filechecker requirement
        
        #corr_array = array(NA, c(matchup$n_levels, matchup$n_prof))
        #corr_array[, matchup$id_prof] = corr
        
        corr_array = ncvar_get(fnc, PARAM_ADJUSTED_NAME)
        
        ### Error
        
        if (PARAM_NAMES[i] == "DOWNWELLING_PAR") {
            ramp_error = 0.05 # 5%
            abs_error = 3e-2 # umol/m²/s # TODO : validate
        } else {
            ramp_error = 0.02 # 2%
            abs_error = 2.5e-5 # W/m²/nm
        }
        
        corr_error_array = pmax(ramp_error*abs(corr_array), abs_error) #TODO : Include applied offset ?
        #corr_error = rep(NA, length(corr))
        #corr_error[which(!is.na(corr))] = -1
        
        #corr_error_array = array(NA, c(matchup$n_levels, matchup$n_prof))
        #corr_error_array[, matchup$id_prof] = corr_error
        
        ### Scientific comment
        
        scientific_comment = str_pad(paste0(PARAM_NAMES[i], " dark correction. Uses JULD to correct drift and SENSOR_TEMP to correct temperature variance. SENSOR_TEMP is reconstructed from the TEMP axis of the core file following [https://doi.org/10.13155/62466]"), 256, "right")
        
        ### Scientific coefficient
        
        #A_total = A_axis_drift_corr[i] + A_axis_corr[i]
        #scientific_coefficient = paste0("A = ", signif(A_total,4),", B = ", signif(B_axis_corr[i],4), ", C = ", signif(C_axis_drift_corr[i],4))
        #if (Q_axis_drift_corr[i] != 0) {
        #    scientific_coefficient = paste0(scientific_coefficient, ", Q = ", signif(Q_axis_drift_corr[i],4))
        #}
        
        ### Scientific equation
        
        #scientific_equation = paste0(PARAM_NAMES[i], "_ADJUSTED = ", PARAM_NAMES[i], " - A - B*SENSOR_TEMP - C*JULD")
        #if (Q_axis_drift_corr[i] != 0) {
        #    scientific_equation = paste0(scientific_equation, " - Q*JULD^2")
        #}
        
        ### comment_dmqc_operator_PRIMARY 
        
        #comment_dmqc_operator_PRIMARY = "PRIMARY | https://orcid.org/0000-0001-9992-5334 | Raphaelle Sauzede, CNRS" 
        
        ### comment_dmqc_operator_PARAM 
        
        #comment_dmqc_operator_PARAM = paste(PARAM_NAMES[i],  "| https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS")
        
        ### HISTORY_SOFTWARE 
        
        #HISTORY_SOFTWARE = "RADM"
        
        ### HISTORY_SOFTWARE_RELEASE
        
        #HISTORY_SOFTWARE_RELEASE = "1.01"
        
        ### write to file
        
        #exit[i] = write_DM(file_out=full_file_name_out, param_name=PARAM_NAMES[i], DATE=DATE, scientific_comment=scientific_comment, 
        #                   scientific_coefficient=scientific_coefficient, scientific_equation=scientific_equation, 
        #                   comment_dmqc_operator_PRIMARY=comment_dmqc_operator_PRIMARY, comment_dmqc_operator_PARAM=comment_dmqc_operator_PARAM, 
        #                   HISTORY_SOFTWARE=HISTORY_SOFTWARE, HISTORY_SOFTWARE_RELEASE=HISTORY_SOFTWARE_RELEASE, param_adjusted=corr_array, 
        #                   param_adjusted_qc=corr_qc_array, param_adjusted_error=corr_error_array)
        
        ncvar_put(fnc, PARAM_ADJUSTED_QC_NAME, corr_qc_array)
        ncvar_put(fnc, PARAM_ADJUSTED_ERROR_NAME, corr_error_array)
       	ncvar_put(fnc, "SCIENTIFIC_CALIB_COMMENT", scientific_comment, start=c(1,id_param,1,id_prof), count=c(256,1,1,1)) 
    }
    
    nc_close(fnc)
    
    return(0)	
    
}

n_cores = detectCores()

res = NULL

#for (WMO in WMO_list) {
for (WMO in c("6901486", "6902880")) {
	print(WMO)
    #files_list = system2("ls", c(paste0("/DATA/correct_RADM/RADM_", WMO, "/RADM_profiles/B*")), stdout = T)
    files_list = system2("ls", c(paste0("/mnt/c/DATA/correct_RADM/RADM_", WMO, "/RADM_profiles/B*")), stdout = T)
    corr_all = mcmapply(corr_file, file_name=files_list, mc.cores=n_cores, SIMPLIFY=FALSE, USE.NAMES=FALSE, MoreArgs=list("WMO"=WMO))
    res[[WMO]] = corr_all
}

#corr_file(file_list[1])

#corr_all = mcmapply(corr_file, file_name=files_list, PROFILE_DATE=date_list, mc.cores=n_cores, SIMPLIFY=FALSE,
#                    MoreArgs=list(path_to_netcdf=path_to_netcdf))

is_ok = rep(NA, length(WMO_list))
for (i in 1:length(is_ok)) {
	is_ok[i] = all(unlist(res[[WMO_list[i]]])==0)
}
