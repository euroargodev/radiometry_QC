require(ncdf4)
require(stringr)
require(stringi)
require(parallel)
require(ggplot2)
require(plyr)
require(viridis)
require(gridExtra)
require(gtools)
require(nortest)

my_menu <- function(choices, title=NULL) { # equivalent to menu() but compatible with non interactive sessions
    
    choice_lines = paste0(as.character(1:length(choices)), 
                          rep(": "),
                          choices)
    cat(title, "", choice_lines, sep="\n")
    
    repeat {
        cat("\nSelection: ")
        
        if (interactive()) {
            input <- as.numeric(readLines(stdin(), 1))
        } else {
            input <- as.numeric(readLines(file("stdin"), 1))
        }

        if ( !is.na(input) ){
            if ( input>=0 & input<=length(choices) ) {
                closeAllConnections()
                return( input )
            }
        }
        
        cat("Choose one of the menu items, or 0 to exit")
    }
}

main_RADM <- function(WMO, index_ifremer, index_greylist, path_to_netcdf, n_cores=detectCores()) {

    cat("Preparing file lists, dates, and greylist...")
    
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
    
    cat("DONE\nImporting and extracting drift data...")
    
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
    
    cat("DONE\nExtracting dark data for method B...")
    
    files_drift_PARALLEL = rep(files_list, 4)
    date_drift_PARALLEL = rep(date_list, 4)
    PARAM_NAMES_drift_PARALLEL = rep(PARAM_NAMES, each=length(files_list))
    
    DRIFT_MATCH = mcmapply(get_profile_match, file_name=files_drift_PARALLEL, param_name=PARAM_NAMES_drift_PARALLEL,
    						path_to_netcdf=path_to_netcdf, PROFILE_DATE=date_drift_PARALLEL, mc.cores=n_cores, USE.NAMES=FALSE, 
    						MoreArgs=list(method="drift_xing"))
    DRIFT_dataf = data.frame("PARAM"=unlist(DRIFT_MATCH[1,]), 
    						"PARAM_Ts"=unlist(DRIFT_MATCH[2,]), 
    					   	"PARAM_date"=unlist(DRIFT_MATCH[3,]), 
    						"PARAM_name"=unlist(DRIFT_MATCH[4,]))
    
    #plot(T_PAR[nna], T_380[nna])
    #lm(T_380[nna] ~ 0 + T_PAR[nna])$coefficients
    
    cat("DONE\nFlagging outliers and greylist...")
    
    ### remove drift outliers
    is_drift_outlier = rep(NA, length(drift_dataf$PARAM))
    is_DRIFT_outlier = rep(NA, length(DRIFT_dataf$PARAM))
    for (param_name in PARAM_NAMES) {
    	param_axis = which(drift_dataf$PARAM_name == param_name)
    
    	quartiles = quantile(drift_dataf$PARAM[param_axis], probs=c(0.25, 0.75))
    	margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
    	outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)
    
    	is_drift_outlier[param_axis] = (drift_dataf$PARAM[param_axis] < outliers_lim[1] | drift_dataf$PARAM[param_axis] > outliers_lim[2])
    
    	### Alternative method
    	param_axis = which(DRIFT_dataf$PARAM_name == param_name)
    
    	quartiles = quantile(DRIFT_dataf$PARAM[param_axis], probs=c(0.25, 0.75))
    	margin = as.numeric(1.5 * (quartiles[2] - quartiles[1]))
    	outliers_lim = c(quartiles[1] - margin, quartiles[2] + margin)
    
    	is_DRIFT_outlier[param_axis] = (DRIFT_dataf$PARAM[param_axis] < outliers_lim[1] | DRIFT_dataf$PARAM[param_axis] > outliers_lim[2])
    }
    
    drift_dataf$is_drift_outlier = is_drift_outlier
    DRIFT_dataf$is_drift_outlier = is_DRIFT_outlier
    
    drift_dataf$PARAM_date_squared = (drift_dataf$PARAM_date)^2
    DRIFT_dataf$PARAM_date_squared = (DRIFT_dataf$PARAM_date)^2
    
    drift_dataf$is_greylisted = mapply(is_greylisted, julian_day=drift_dataf$PARAM_date,
    								PARAMETER_NAME=as.character(drift_dataf$PARAM_name),
    								MoreArgs=list(WMO=WMO))
    DRIFT_dataf$is_greylisted = mapply(is_greylisted, julian_day=DRIFT_dataf$PARAM_date,
                                       PARAMETER_NAME=as.character(DRIFT_dataf$PARAM_name),
                                       MoreArgs=list(WMO=WMO))
    
    cat("DONE\n")
    
    do_quadratic_fit = c(F, F, F, F)

    repeat { # loop until the user finds an adequate fit or quits
        
        cat("Fitting both methods...")
        
        A_axis_drift = rep(NA, 4)
        B_axis_drift = rep(NA, 4)
        C_axis_drift = rep(NA, 4)
        Q_axis_drift = rep(0, 4)
        drift_dataf_5C = drift_dataf
        fitted_coeff_drift = NULL
        A_axis_DRIFT = rep(NA, 4)
        B_axis_DRIFT = rep(NA, 4)
        C_axis_DRIFT = rep(NA, 4)
        DRIFT_dataf_5C = DRIFT_dataf
        fitted_coeff_DRIFT = NULL
        Q_axis_DRIFT = rep(0, 4)
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
        }
        for (i in 1:4) {
        	### Alternative drift
        	subset_PAR = which(DRIFT_dataf$PARAM_name == PARAM_NAMES[i])
        	subset_fit = which(DRIFT_dataf$PARAM_name == PARAM_NAMES[i] & !DRIFT_dataf$is_greylisted
        	                   & !DRIFT_dataf$is_drift_outlier)
        	if (length(subset_fit) == 0) { next }	
        	
        	if (do_quadratic_fit[i]) {
        	    fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date + PARAM_date_squared, data=DRIFT_dataf, subset=subset_fit)
        	    Q_axis_DRIFT[i] = fit_ABC$coefficients[4]	
        	} else {
        	    fit_ABC = lm(PARAM ~ PARAM_Ts + PARAM_date, data=DRIFT_dataf, subset=subset_fit)
        	}
        
        	fitted_coeff_DRIFT[[PARAM_NAMES[i]]] = fit_ABC$coefficients
        	A_axis_DRIFT[i] = fit_ABC$coefficients[1]	
        	B_axis_DRIFT[i] = fit_ABC$coefficients[2]	
        	C_axis_DRIFT[i] = fit_ABC$coefficients[3]	
        
        	DRIFT_dataf_5C$PARAM[subset_PAR] = DRIFT_dataf_5C$PARAM[subset_PAR] - B_axis_DRIFT[i] * (DRIFT_dataf_5C$PARAM_Ts[subset_PAR] - 5)
        }
        
        cat("DONE\nCreating drift plots...")
        
        good_drift = which(!drift_dataf$is_greylisted & !drift_dataf$is_drift_outlier)
        range_time = seq(min(c(drift_dataf$PARAM_date[good_drift], date_list)), max(c(drift_dataf$PARAM_date[good_drift], date_list)), length.out=100 )
        
        data_fit_drift = data.frame(
        	PARAM_name = rep(PARAM_NAMES, each=2),
        	x = rep(range_time, 4),
        	y = rep(A_axis_drift, each=2) + rep(B_axis_drift, each=2) * 5
        		+ rep(C_axis_drift, each=2) * rep(range_time, 4) 
        		+ rep(Q_axis_drift, each=2) * rep(range_time^2, 4) 
        )
        data_fit_DRIFT = data.frame(
        	PARAM_name = rep(PARAM_NAMES, each=2),
        	x = rep(range_time, 4),
        	y = rep(A_axis_DRIFT, each=2) + rep(B_axis_DRIFT, each=2) * 5 
        	    + rep(C_axis_DRIFT, each=2) * rep(range_time, 4) 
        	    + rep(Q_axis_DRIFT, each=2) * rep(range_time^2, 4) 
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
        g6 = ggplot(na.omit(DRIFT_dataf), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
        	geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
            geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
        	scale_color_viridis() +
        	facet_wrap(~PARAM_name, scale="free_y") +
            labs (x="Julian day", y="Irradiance", colour="Sensor\ntemperature", 
                  title="Dark data extracted from profiles")
        g7 = ggplot(na.omit(DRIFT_dataf_5C), aes(x=PARAM_date, y=PARAM, color=PARAM_Ts)) +
            geom_point(data=function(x){x[!x$is_greylisted & !x$is_drift_outlier, ]}) +
            geom_point(data=function(x){x[x$is_greylisted & !x$is_drift_outlier, ]}, color="red") +
        	scale_color_viridis() +
        	facet_wrap(~PARAM_name, scale="free_y") +
            labs (x="Julian day", y="Irradiance", colour="Sensor\ntemperature", 
                  title="Dark data extracted from profiles,\nfitted to time and temperature and adjusted to 5°C")
        	
        #g4_fit = g4 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
        g5_fit = g5 + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
        g5_fit_2 = g5_fit + geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")
        
        #g6_fit = g6 + geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")
        g7_fit = g7 + geom_line(data=data_fit_DRIFT, mapping=aes(x=x,y=y), color="black")
        g7_fit_2 = g7_fit + geom_line(data=data_fit_drift, mapping=aes(x=x,y=y), color="red")
        
        if ( dim(drift_dataf)[1] != 0 ) {
            x11(xpos=0, ypos=0)
            plot(g4)
            x11(xpos=0, ypos=100)
            plot(g5_fit)
            x11(xpos=0, ypos=200)
            plot(g5_fit_2)
            
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
        if ( dim(DRIFT_dataf)[1] != 0 ) {
            x11(xpos=-1, ypos=0)
            plot(g6)
            x11(xpos=-1, ypos=100)
            plot(g7_fit)
            x11(xpos=-1, ypos=200)
            plot(g7_fit_2)
            
            plot_name=paste0(WMO, "_regr_drift_B.png")
            png(filename=plot_name, width=600, height=600)
            plot(g7_fit)
            dev.off()
            
            plot_name=paste0(WMO, "_regr_drift_uncorrected_B.png")
            png(filename=plot_name, width=600, height=600)
            plot(g6)
            dev.off()
        } else {
            cat("The alternative drift method found no valid dark data points...")
        }
        
        cat("DONE\n")
        
        choice = my_menu(title = "What is the next step ? (0 to abandon and quit)",
                         choices = c("Continue with correction from method A",
                                     "Continue with correction from method B",
                                     "Change to quadratic fits for some parameters",
                                     "Continue without a drift correction"))
        switch (choice + 1,
            return(0),
            
            {fitted_coeff_drift_corr = fitted_coeff_drift
            A_axis_drift_corr = A_axis_drift
            B_axis_drift_corr = B_axis_drift
            C_axis_drift_corr = C_axis_drift
            Q_axis_drift_corr = Q_axis_drift
            break},
            
            {fitted_coeff_drift_corr = fitted_coeff_DRIFT
            A_axis_drift_corr = A_axis_DRIFT
            B_axis_drift_corr = B_axis_DRIFT
            C_axis_drift_corr = C_axis_DRIFT
            Q_axis_drift_corr = Q_axis_DRIFT
            break},
            
            { },
            
            {fitted_coeff_drift_corr = NULL
            A_axis_drift_corr = c(0, 0, 0, 0)
            B_axis_drift_corr = c(0, 0, 0, 0)
            C_axis_drift_corr = c(0, 0, 0, 0)
            Q_axis_drift_corr = c(0, 0, 0, 0)
            break}
        )
        
        do_quadratic_fit = c(F, F, F, F)
        for (i in 1:4) {
            choice = my_menu(title = paste("Use a quadratic fit for", PARAM_NAMES[i], "?"),
                          choices = c("Yes", "No"))
            do_quadratic_fit[i] = (choice == 1)
        }
        
        while(dev.cur() > 1) {dev.off()}
    }
    
    while(dev.cur() > 1) {dev.off()}
    
    
    
    ###############################################
    #### Start of radiometry/Ts matches extraction
    ###############################################
    
    cat("Extracting day and night profiles...")
    
    files_night_PARALLEL = rep(files_night, 4)
    date_night_PARALLEL = rep(date_night, 4)
    PARAM_NAMES_night_PARALLEL = rep(PARAM_NAMES, each=length(files_night))
    A_axis_drift_night_PARALLEL = rep(A_axis_drift_corr, each=length(files_night))
    C_axis_drift_night_PARALLEL = rep(C_axis_drift_corr, each=length(files_night))
    Q_axis_drift_night_PARALLEL = rep(Q_axis_drift_corr, each=length(files_night))
    
    files_day_PARALLEL = rep(files_day, 4)
    date_day_PARALLEL = rep(date_day, 4)
    PARAM_NAMES_day_PARALLEL = rep(PARAM_NAMES, each=length(files_day))
    A_axis_drift_day_PARALLEL = rep(A_axis_drift_corr, each=length(files_day))
    C_axis_drift_day_PARALLEL = rep(C_axis_drift_corr, each=length(files_day))
    Q_axis_drift_day_PARALLEL = rep(Q_axis_drift_corr, each=length(files_day))
    
    NIGHT_MATCH = mcmapply(get_profile_match, file_name=files_night_PARALLEL, param_name=PARAM_NAMES_night_PARALLEL,
    						PROFILE_DATE=date_night_PARALLEL, path_to_netcdf=path_to_netcdf, drift_A=A_axis_drift_night_PARALLEL,
    						drift_C=C_axis_drift_night_PARALLEL, drift_Q=Q_axis_drift_night_PARALLEL,
    						mc.cores=n_cores, USE.NAMES=FALSE)
    DAY_MATCH = mcmapply(get_profile_match, file_name=files_day_PARALLEL, param_name=PARAM_NAMES_day_PARALLEL,
    						PROFILE_DATE=date_day_PARALLEL, path_to_netcdf=path_to_netcdf, drift_A=A_axis_drift_day_PARALLEL,
    						drift_C=C_axis_drift_day_PARALLEL, drift_Q=Q_axis_drift_day_PARALLEL,
    						mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(method="day"))
    
    cat("DONE\nFlagging greylisted profiles and outliers in day profiles...")
    
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
    					   "PARAM_date"=unlist(NIGHT_MATCH[3,]), "PARAM_name"=unlist(NIGHT_MATCH[4,]),
    					   "PARAM_pres"=unlist(NIGHT_MATCH[5,]))
    DAY_dataf = data.frame("PARAM"=unlist(DAY_MATCH[1,]), 
    						"PARAM_Ts"=unlist(DAY_MATCH[2,]), 
    					   	"PARAM_date"=unlist(DAY_MATCH[3,]), 
    						"PARAM_name"=unlist(DAY_MATCH[4,]),
    						"is_dark_outlier"=is_dark_outlier_full,
    						"PARAM_pres"=unlist(DAY_MATCH[7,]))
    
    PAR_dataf$is_greylisted = mapply(is_greylisted, julian_day=PAR_dataf$PARAM_date,
    								PARAMETER_NAME=as.character(PAR_dataf$PARAM_name),
    								MoreArgs=list(WMO=WMO))
    DAY_dataf$is_greylisted = mapply(is_greylisted, julian_day=DAY_dataf$PARAM_date,
    								PARAMETER_NAME=as.character(DAY_dataf$PARAM_name),
    								MoreArgs=list(WMO=WMO))
    
    
    cat("DONE\nExtracting full Ts range...")
    
    TEST_MATCH = mcmapply(get_profile_match, file_name=files_drift_PARALLEL, param_name=PARAM_NAMES_drift_PARALLEL,
                          path_to_netcdf=path_to_netcdf, PROFILE_DATE=date_drift_PARALLEL, mc.cores=n_cores, USE.NAMES=FALSE, 
                          MoreArgs=list(method="test"))
    TEST_dataf = data.frame("PARAM"=unlist(TEST_MATCH[1,]), 
                            "PARAM_Ts"=unlist(TEST_MATCH[2,]), 
                            "PARAM_date"=unlist(TEST_MATCH[3,]), 
                            "PARAM_name"=unlist(TEST_MATCH[4,]))
    
    #T_PAR = TEST_dataf$PARAM[which(TEST_dataf$PARAM_name=="DOWNWELLING_PAR")]
    #T_380 = TEST_dataf$PARAM[which(TEST_dataf$PARAM_name=="DOWN_IRRADIANCE380")]
    #T_412 = TEST_dataf$PARAM[which(TEST_dataf$PARAM_name=="DOWN_IRRADIANCE412")]
    #T_490 = TEST_dataf$PARAM[which(TEST_dataf$PARAM_name=="DOWN_IRRADIANCE490")]
    
    #nna = which(!is.na(T_PAR) & !is.na(T_380) & !is.na(T_412) & !is.na(T_490))
    
    Ts_range = range(c(TEST_dataf$PARAM_Ts), na.rm=T)
    
    cat("DONE\n")
    
    ###########################
    ### Start of Ts fitting
    ###########################
    
    pres_cutoff = -Inf
    
    repeat {
    
        cat("Fitting Ts to radiometry parameters...")
        
        fitted_coeff = NULL
        A_axis = rep(NA, 4)
        B_axis = rep(NA, 4)
        fitted_coeff_day = NULL
        A_axis_day = rep(NA, 4)
        B_axis_day = rep(NA, 4)
        for (i in 1:length(PARAM_NAMES)) {
            
            ## night
        	subset_PAR = which(PAR_dataf$PARAM_name==PARAM_NAMES[i] & !is.na(PAR_dataf$PARAM_Ts)
        						& !PAR_dataf$is_greylisted & PAR_dataf$PARAM_pres>pres_cutoff)
        
        	fit_AB = lm(PARAM ~ PARAM_Ts, data=PAR_dataf, subset=subset_PAR) 
        	
        	fitted_coeff[[PARAM_NAMES[i]]] = fit_AB$coefficients
        	A_axis[i] = fit_AB$coefficients[[1]]
        	B_axis[i] = fit_AB$coefficients[[2]]
        	
        	## day
        	subset_DAY = which(DAY_dataf$PARAM_name==PARAM_NAMES[i] & !is.na(DAY_dataf$PARAM_Ts)
        						& !DAY_dataf$is_greylisted & !DAY_dataf$is_dark_outlier
        						& DAY_dataf$PARAM_pres>pres_cutoff)
        
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
        
        
        
        cat("DONE\nCreating Ts/Irr plots...")
    
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
        
        g1 = ggplot(na.omit(PAR_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_pres, group=PARAM_name)) +
            geom_point(data=function(x){x[!x$is_greylisted & x$PARAM_pres>pres_cutoff, ]}) +
            #geom_point(data=function(x){x[x$is_greylisted, ]}, color="red") +
            scale_color_viridis(trans='reverse') +
            facet_wrap(~PARAM_name, scale="free_y") +
            labs (x="Sensor temperature", y="Irradiance", colour="Pressure (db)", 
                  title="Data from night profiles")
        g1_day = ggplot(na.omit(DAY_dataf), aes(x=PARAM_Ts, y=PARAM, color=PARAM_pres, group=PARAM_name)) +
            geom_point(data=function(x){x[!x$is_greylisted & !x$is_dark_outlier & x$PARAM_pres>pres_cutoff, ]}) +
            #geom_point(data=function(x){x[x$is_greylisted & !x$is_dark_outlier, ]}, color="red") +
            scale_color_viridis(trans='reverse') +
            facet_wrap(~PARAM_name, scale="free_y") +
            labs (x="Sensor temperature", y="Irradiance", colour="Pressure (db)", 
                  title="Data selected from day profiles")
        
        g2 = g1 + geom_line(data=data_fit, mapping=aes(x=x,y=y), color="red")
        g3 = g2 + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="black")
        
        g2_day = g1_day + geom_line(data=data_fit_day, mapping=aes(x=x,y=y), color="black")
        g3_day = g2_day + geom_line(data=data_fit, mapping=aes(x=x,y=y), color="red") 
        
        x11(xpos=0, ypos=0)
        plot(g2)
        x11(xpos=0, ypos=-1)
        plot(g3)
        x11(xpos=-1, ypos=0)
        plot(g2_day)
        x11(xpos=-1, ypos=-1)
        plot(g3_day)
        
        plot_name=paste0(WMO, "_regr_day.png")
        png(filename=plot_name, width=600, height=600)
        plot(g3_day)
        dev.off()
        
        plot_name=paste0(WMO, "_regr_night.png")
        png(filename=plot_name, width=600, height=600)
        plot(g3)
        dev.off()
        
        cat("DONE\n")
        
        choice = my_menu(title = "Which correction should be used ? (0 to abandon and quit)",
                         choices = c("Continue with correction from night method",
                                     "Continue with correction from day method",
                                     "Add or change the pressure cutoff for data selection"))
        
        switch (choice + 1,
                return(0),
                
                {fitted_coeff_corr = fitted_coeff
                A_axis_corr = A_axis
                B_axis_corr = B_axis
                break},
                
                {fitted_coeff_corr = fitted_coeff_day
                A_axis_corr = A_axis_day
                B_axis_corr = B_axis_day
                break},
                
                { }
        )
        
        cat("New pressure cutoff in decibars: ")
        
        if (interactive()) {
            pres_cutoff <- as.numeric(readLines(stdin(), 1))
        } else {
            pres_cutoff <- as.numeric(readLines(file("stdin"), 1))
        }
        
        while(dev.cur() > 1) {dev.off()}
 
    }
        
    while(dev.cur() > 1) {dev.off()}
    
    choice = my_menu(title = "What QC flag is the DM allowed to have at best in the dark section of profiles ? (0 to abandon and quit)",
                     choices = c("Good",
                                 "Probably good",
                                 "Probably bad",
                                 "Bad"))
    if (choice == 0) { return(0) }
    
    best_QC = choice
    
    #### Deal with missing dates
    
    prof_num_date_miss = substr(prof_id[subset_date_miss], 11, 13)
    if ( length(prof_num_date_miss) != 0 ) {
        
        question = paste0("Profiles with missing date : ", paste(prof_num_date_miss, collapse = ", "), 
                          "\nWhat should be done with them ? (0 to abandon and quit)")
        
        choice = my_menu(title = question,
                         choices = c("Attempt to find a suited date for DM correction",
                                     "Ignore the profiles (just copy them)"))
        
        date_list_missing = rep(NA, length(prof_num_date_miss))
        
        if ( choice == 0 ) {
            return(0)
        } else if (choice == 1) {
            for (i in 1:length(prof_num_date_miss)) {
                
                date_prev = date_list[which(as.numeric(substr(prof_id[subset], 11, 13)) == (as.numeric(prof_num_date_miss[1]) - 1) )]
                date_next = date_list[which(as.numeric(substr(prof_id[subset], 11, 13)) == (as.numeric(prof_num_date_miss[1]) + 1) )]
                
                question2 = c("For profile ", prof_num_date_miss[i], " :\n    Date of preceding profile : ", NA, "\n    Date of following profile : ", NA, 
                              "\nWhat should be done with this profile ? (If any date is unknown, do not use it)")
                
                if (length(date_prev) != 0) {
                    question2[4] = as.character(as.Date(date_prev, origin=as.Date("1950-01-01", tz="UTC")))
                } else {
                    question2[4] = "Unknown"
                }
                if (length(date_next) != 0) {
                    question2[6] = as.character(as.Date(date_next, origin=as.Date("1950-01-01", tz="UTC")))
                } else {
                    question2[6] = "Unknown"
                }
                
                question2 = paste(question2, collapse="")
                choice2 = my_menu(title = question2,
                                  choices = c("Use the preceding date",
                                              "Use the following date",
                                              "Use the interpolation of the dates",
                                              "Ignore this profile (just copy it)"))
                switch (choice2 + 1,
                        return(0),
                        { date_list_missing[i] = date_prev },
                        { date_list_missing[i] = date_next },
                        { date_list_missing[i] = (date_prev + date_next) / 2 }
                )
            }
        }
        
        # move the profiles with no alternative date to the just_copy list
        subset_copy = c(subset_copy, subset_date_miss[ which(is.na(date_list_missing)) ])
        subset_date_miss = subset_date_miss[ which(!is.na(date_list_missing)) ]
        date_list_missing = date_list_missing[ which(!is.na(date_list_missing)) ]
        
    } else {
        cat("No profiles with missing date\n")
    }
    
    ### Mention missing positions to user
    
    if ( length(subset_pos_miss) != 0 ) {
        cat(paste0("Profiles with missing position : ", paste(substr(prof_id[subset_pos_miss], 11, 13), collapse = ", "), "\n"))
    } else {
        cat("No profiles with missing position\n")
    }
    
    
    
    ####################################
    ### Apply correction to netcdf
    ####################################
    
    date_update = Sys.time()
    DATE = stri_datetime_format(date_update, format="uuuuMMddHHmmss", tz="UTC")
    
    corr_file <- function(file_name, PROFILE_DATE, path_to_netcdf, just_copy=FALSE) {
    
    	path_sep = unlist(strsplit(file_name, "/"))
    	
    	if (!just_copy) { str_sub(path_sep[4], 2, 2) = "D" }

    	file_name_out = paste(path_sep[1], path_sep[2], path_sep[3], "RADM/RADM_profiles", path_sep[4], sep="/")

    	full_file_name_out = paste(path_to_netcdf, file_name_out, sep="")
    	full_file_name_in = paste(path_to_netcdf, file_name, sep="")
    	
    	system2("cp", c(full_file_name_in, full_file_name_out))
    	
    	if (just_copy) { return(rep(0, 4)) }
        
    	exit = rep(NA, 4)
    	
    	for (i in 1:length(PARAM_NAMES)) {
    	    
    	    matchup = get_Ts_match(file_name=file_name, path_to_netcdf=path_to_netcdf, PARAM_NAME=PARAM_NAMES[i])
    	    
    	    match_day = get_profile_match(file_name=file_name, param_name=PARAM_NAMES[i], path_to_netcdf=path_to_netcdf, 
    	                                  PROFILE_DATE=PROFILE_DATE, method="day", drift_A=A_axis_drift_corr[i], 
    	                                  drift_C=C_axis_drift_corr[i], drift_Q=Q_axis_drift_corr[i])
    	    
    	    ### Flags
    	    
    	    corr_qc = matchup$PARAM_QC
    	    
    	    corr_qc[match_day$MATCH_which] = rep(as.character(best_QC), length(match_day$MATCH_which)) #if which==NULL --> no change
    	    corr_qc[which(matchup$PARAM_QC=="3" | matchup$PARAM_QC=="4")] = "4"
    	    corr_qc[which(matchup$PRES_B_QC=="3" | matchup$PRES_B_QC=="4")] = "4"
    	    corr_qc[which(is.na(matchup$Ts) & !is.na(matchup$PARAM))] = "4"
    	    corr_qc_word = paste(corr_qc, collapse="")	
    	    
    	    blank_qc = paste(rep(" ", matchup$n_levels), collapse="")
    	    corr_qc_array = rep(blank_qc, matchup$n_prof)
    	    corr_qc_array[matchup$id_prof] = corr_qc_word
    		
    	    ### Adjusted parameter
    	    
    		param_undrifted = as.numeric( matchup$PARAM - A_axis_drift_corr[i] - C_axis_drift_corr[i] * PROFILE_DATE
    										- Q_axis_drift_corr[i] * PROFILE_DATE^2 )
    		
    		corr = param_undrifted - ( A_axis_corr[i] + B_axis_corr[i] * matchup$Ts )
    		corr[which(corr_qc == "4")] = NA #filechecker requirement
    		
    		corr_array = array(NA, c(matchup$n_levels, matchup$n_prof))
    		corr_array[, matchup$id_prof] = corr
    		
    		### Error
    		
    		#corr_error = 0.2*abs(corr) #TODO
    		corr_error = rep(NA, length(corr))
    		corr_error[which(!is.na(corr))] = -1
    		
    		corr_error_array = array(NA, c(matchup$n_levels, matchup$n_prof))
    		corr_error_array[, matchup$id_prof] = corr_error
            
    		### Scientific comment
            
            scientific_comment = paste0(PARAM_NAMES[i], " dark correction. Uses JULD to correct drift and SENSOR_TEMP to correct temperature variance. SENSOR_TEMP is reconstructed from the TEMP axis of the core file following [https://doi.org/10.1117/12.2504241] with delta_t=60s and k=12h^-1")
    		### Scientific coefficient
            
            A_total = A_axis_drift_corr[i] + A_axis_corr[i]
            scientific_coefficient = paste0("A = ", signif(A_total,4),", B = ", signif(B_axis_corr[i],4), ", C = ", signif(C_axis_drift_corr[i],4))
            if (Q_axis_drift_corr[i] != 0) {
                scientific_coefficient = paste0(scientific_coefficient, ", Q = ", signif(Q_axis_drift_corr[i],4))
            }
    		
    		### Scientific equation
            
            scientific_equation = paste0(PARAM_NAMES[i], "_ADJUSTED = ", PARAM_NAMES[i], " - A - B*SENSOR_TEMP - C*JULD")
            if (Q_axis_drift_corr[i] != 0) {
                scientific_equation = paste0(scientific_equation, " - Q*JULD^2")
            }
    		
    		### comment_dmqc_operator_PRIMARY 
            
            comment_dmqc_operator_PRIMARY = "PRIMARY | https://orcid.org/0000-0001-9992-5334 | Raphaelle Sauzede, CNRS" 
            
    		### comment_dmqc_operator_PARAM 
            
            comment_dmqc_operator_PARAM = paste(PARAM_NAMES[i],  "| https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS")
    		
    		### HISTORY_SOFTWARE 
            
            HISTORY_SOFTWARE = "RADM"
    		
    		### HISTORY_SOFTWARE_RELEASE
    		
            HISTORY_SOFTWARE_RELEASE = "1.00"
    		
    		### write to file
    		
    		exit[i] = write_DM(file_out=full_file_name_out, param_name=PARAM_NAMES[i], DATE=DATE, scientific_comment=scientific_comment, 
    		                   scientific_coefficient=scientific_coefficient, scientific_equation=scientific_equation, 
    		                   comment_dmqc_operator_PRIMARY=comment_dmqc_operator_PRIMARY, comment_dmqc_operator_PARAM=comment_dmqc_operator_PARAM, 
    		                   HISTORY_SOFTWARE=HISTORY_SOFTWARE, HISTORY_SOFTWARE_RELEASE=HISTORY_SOFTWARE_RELEASE, param_adjusted=corr_array, 
    		                   param_adjusted_qc=corr_qc_array, param_adjusted_error=corr_error_array)

    	}
    	
    	return(exit)	
    		
    }
    
    cat("Applying correction to netcdf files...")
    
    #corr_file(files_list[1], PROFILE_DATE=date_list[1], path_to_netcdf=path_to_netcdf)
    
    corr_all = mcmapply(corr_file, file_name=files_list, PROFILE_DATE=date_list, mc.cores=n_cores, SIMPLIFY=FALSE,
    							MoreArgs=list(path_to_netcdf=path_to_netcdf))
    
    if ( length(subset_date_miss) != 0 ) {
        corr_date_miss = mcmapply(corr_file, file_name=files[subset_date_miss], PROFILE_DATE=date_list_missing, mc.cores=n_cores, SIMPLIFY=FALSE,
                                    MoreArgs=list(path_to_netcdf=path_to_netcdf))
    }
    
    corr_copy = mcmapply(corr_file, file_name=files[subset_copy], mc.cores=n_cores, SIMPLIFY=FALSE,
                         MoreArgs=list(path_to_netcdf=path_to_netcdf, PROFILE_DATE=NA, just_copy=TRUE))
    

    cat("DONE\n")
    
    return(0)
}


