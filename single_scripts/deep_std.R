###################################################################
### compute radiometry std at depth for error estimation
###################################################################

### import pathways
source("~/Documents/radiometry/pathways.R")
### Source the DM functions and subfunctions, and libraries
suppressPackageStartupMessages({
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
    library(tidyr)
    library(dplyr)
})

source(paste0(path_to_source, "possol.R"))
source(paste0(path_to_source, "sensor_temp.R"))
source(paste0(path_to_source, "main_RADM.R"))
source(paste0(path_to_source, "get_matches.R"))
source(paste0(path_to_source, "write_DM.R"))
source(paste0(path_to_source, "increment_N_CALIB.R"))
source(paste0(path_to_source, "RT_QC_radiometry_function_oao_2.R"))
source(paste0(path_to_source, "plots.R"))

### import tables
index_ifremer = read.table(path_to_index_ifremer, sep=",", header = T)
index_greylist = read.csv(path_to_index_greylist, sep = ",")

get_deep_sd <- function(WMO, n_cores=detectCores()) {
    
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
    #WMO = "6901473" # again drift data missing, maybe other sensor issues
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
    
    subset = which( wod==WMO & substr(prof_id,14,14)!="D" & !is.na(prof_date) & grepl("DOWNWELLING_PAR", variables) )
    subset_date_miss = which( wod==WMO & substr(prof_id,14,14)!="D" & is.na(prof_date) & grepl("DOWNWELLING_PAR", variables) )
    subset_copy = which( wod==WMO & ( substr(prof_id,14,14)=="D" | !grepl("DOWNWELLING_PAR", variables) ) )
    
    subset_pos_miss = which( wod==WMO & substr(prof_id,14,14)!="D" & grepl("DOWNWELLING_PAR", variables) & (is.na(lat) | is.na(lon)) )
    
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
    
    
    DRIFT = numeric(0)
    DRIFT_Ts = numeric(0)
    DRIFT_date = numeric(0)
    DRIFT_name = character(0)
    for (param_name in PARAM_NAMES) {
        irr = ncvar_get(fnc_B, param_name)
        
        drift_B = which(!is.na(irr) & code==290) #290 may not be the only correct code
        
        if ( length(drift_B)==0 ) {
            next
        }
        
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
    
    ### remove drift outliers
    is_drift_outlier = rep(NA, length(drift_dataf$PARAM))
    for (param_name in PARAM_NAMES) {
        param_axis = which(drift_dataf$PARAM_name == param_name)
        
        quartiles = quantile(drift_dataf$PARAM[param_axis], probs=c(0.25, 0.75))
        margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
        outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)
        
        is_drift_outlier[param_axis] = (drift_dataf$PARAM[param_axis] < outliers_lim[1] | drift_dataf$PARAM[param_axis] > outliers_lim[2])
    }
    
    drift_dataf$is_drift_outlier = is_drift_outlier
    
    drift_dataf$PARAM_date_squared = (drift_dataf$PARAM_date)^2
    
    drift_dataf$is_greylisted = mcmapply(is_greylisted, julian_day=drift_dataf$PARAM_date,
                                         PARAMETER_NAME=as.character(drift_dataf$PARAM_name),
                                         MoreArgs=list(WMO=WMO), mc.cores=n_cores)
    
    do_quadratic_fit = c(F, F, F, F)

    A_axis_drift = rep(NA, 4)
    B_axis_drift = rep(NA, 4)
    C_axis_drift = rep(NA, 4)
    Q_axis_drift = rep(0, 4)
    drift_dataf_5C = drift_dataf
    drift_dataf_corr = drift_dataf
    fitted_coeff_drift = NULL
    remain_sd = rep(NA, 4)
    for (i in 1:4) {
        subset_PAR = which(drift_dataf$PARAM_name == PARAM_NAMES[i])
        subset_fit = which(drift_dataf$PARAM_name == PARAM_NAMES[i] & !drift_dataf$is_greylisted
                           & !drift_dataf$is_drift_outlier)
        if ( length(subset_fit)==0 ) { next }
        
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
        drift_dataf_corr$PARAM[subset_PAR] = drift_dataf_5C$PARAM[subset_PAR] - C_axis_drift[i] * drift_dataf_corr$PARAM_date[subset_PAR]
        
        remain_sd[i] = sd(drift_dataf_corr$PARAM[subset_fit])
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
        
    g4 = ggplot(na.omit(drift_dataf), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
        geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
        geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
        scale_color_viridis() +
        facet_wrap(~PARAM_name, scale="free_y") +
        labs (x="Julian day", y="Irradiance", colour="Closest drift\ntemperature", 
                  title="Irradiance measurements during float drift")
    g5 = ggplot(na.omit(drift_dataf_5C), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
        geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
        geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
        #geom_point(data=function(x){x[x$is_drift_outlier & !x$is_greylisted, ]}, color="red", shape=4) +
        scale_color_viridis() +
        facet_wrap(~PARAM_name, scale="free_y") +
        labs (x="Julian day", y="Irradiance", colour="Closest drift\ntemperature",
              title="Irradiance measurements during float drift,\nfitted to time and temperature and adjusted to 5°C")
        
        
    #g4_fit = g4 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
    g5_fit = g5 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
        
    if ( dim(drift_dataf)[1] != 0 ) {
        #x11(xpos=0, ypos=0)
        #plot(g4)
        #x11(xpos=0, ypos=100)
        #plot(g5_fit)
        #x11(xpos=0, ypos=200)
        #plot(g5_fit_2)
        
        plot_name=paste0(WMO, "_regr_drift_A.png")
        png(filename=plot_name, width=600, height=600)
        plot(g5_fit)
        dev.off()
        
        plot_name=paste0(WMO, "_regr_drift_uncorrected_A.png")
        png(filename=plot_name, width=600, height=600)
        plot(g4)
        dev.off()
    } else {
        cat("No valid drift data points were found...")
    }

    return(remain_sd)
}

WMO_list = c("6901655", "6902735", "6902826", "6902954", "6901004", "6901480", "6901486", "6901764", "6901862", "6901527",
             "6901647", "6902547", "6901474", "6901489", "6901654", "6901657", "6901766", "6901770", "6901771", "6901865",
             "6903197", "6901580", "6901862", "6902732", "7900591", "6901492", "6901519", "6901573", "6901575", "6901767",
             "6901768", "6901773", "6902738", "6902828")

dataf_sd = data.frame("WMO"=WMO_list, "PAR_sd"=numeric(length(WMO_list)), "I380_sd"=numeric(length(WMO_list)), "I412_sd"=numeric(length(WMO_list)), "I490_sd"=numeric(length(WMO_list)))

for (i in 1:length(WMO_list)) {
    dataf_sd[i, 2:5] = get_deep_sd(WMO_list[i])
}

hist(dataf_sd$PAR_sd, xlab="PAR standard deviation (umol/m²/s)", main="Histogram of the remaining SD in drift data after a\n multi-linear regression (34 hand-selected floats)")
hist(dataf_sd$I380_sd, xlab="Ed(380) standard deviation (W/m²/nm)", main="Histogram of the remaining SD in drift data after a\n multi-linear regression (34 hand-selected floats)")
hist(dataf_sd$I412_sd, xlab="Ed(412) standard deviation (W/m²/nm)", main="Histogram of the remaining SD in drift data after a\n multi-linear regression (34 hand-selected floats)")
hist(dataf_sd$I490_sd, xlab="Ed(490) standard deviation (W/m²/nm)", main="Histogram of the remaining SD in drift data after a\n multi-linear regression (34 hand-selected floats)")

WMO_list[which(dataf_sd$PAR_sd>0.01)]
hist(dataf_sd$PAR_sd[which(dataf_sd$PAR_sd<0.01)], xlab="PAR standard deviation (umol/m²/s)", main="Histogram of the remaining SD in drift data after a\n multi-linear regression (30 hand-selected floats)")
