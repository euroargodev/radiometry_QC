sensor_temp <- function(TEMP, PRES_TEMP, PRES_PARAM, material="PEEK") {
	
    if (material=="PEEK") {
	    #k = 0.19/60 # s^-1
	    #delta_t = 54 # s
	    k = 12/3600 # s^-1
	    delta_t = 60 # s
    } else if (material == "Aluminium") {
        k = 0.44/60 # s^-1
        delta_t = 15 # s
    } else {
        return(rep(NA, length(PRES_PARAM)))
    }
	
	asc_speed = 0.1 # m/s

	Tw = rev(TEMP[!is.na(TEMP) & !is.na(PRES_TEMP)])
	Tw_pres = rev(PRES_TEMP[!is.na(TEMP) & !is.na(PRES_TEMP)]) #pressure associated to Tw measurements
	Tw_time = (max(Tw_pres) - Tw_pres) / asc_speed #assuming 1m/db for now, time associated to Tw measurements

	n_time = length(Tw)
	
	if (n_time < 2) {
		return(rep(NA, length(PRES_PARAM)))
	}
	
	Ts = rep(NA, n_time)
	Ts[1] = Tw[1]
	Ts_time = Tw_time + delta_t # time associated to Ts estimates
	Ts_pres = max(Tw_pres) - Ts_time*asc_speed # pressure associated to Ts estimates

	for (n in 2:n_time) {
		Ts[n] = Ts[n-1] + k * (Tw_time[n]-Tw_time[n-1]) * (Tw[n-1]-Ts[n-1])
	}
	
	fitted_Ts = approx(Ts_pres, Ts, xout=PRES_PARAM)$y
	
	return(fitted_Ts)
}
