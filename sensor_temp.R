sensor_temp <- function(TEMP, PRES_TEMP, PRES_PARAM, material="PEEK", MTIME_PARAM=NULL) {
  
  if (material=="PEEK") {
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
  
  n_time = length(Tw)
  
  if (n_time < 2) {
    return(rep(NA, length(PRES_PARAM)))
  }
  
  #pressure associated to Tw measurements
  Tw_pres = rev(PRES_TEMP[!is.na(TEMP) & !is.na(PRES_TEMP)])
  
  if ( is.null(MTIME_PARAM) ) {
    #assuming 1m/db, time associated to Tw measurements
    Tw_time = (max(Tw_pres) - Tw_pres) / asc_speed
  } else {
    Tw_time = approx(rev(PRES_PARAM), rev(MTIME_PARAM)*86400, xout=Tw_pres)$y
    
    # If necessary, extrapolate the remaining of the axis with 0.1 m/s assumption
    interp_range = range(which(!is.na(Tw_time)))
    if ( interp_range[1] != 1 ) {
      Tw_time[1:(interp_range[1] - 1)] = Tw_time[interp_range[1]] + 
        (Tw_pres[interp_range[1]] - Tw_pres[1:(interp_range[1] - 1)]) / asc_speed
    }
    if ( interp_range[2] != length(Tw_time) ) {
      Tw_time[(interp_range[2] + 1):length(Tw_time)] = Tw_time[interp_range[2]] + 
        (Tw_pres[interp_range[2]] - Tw_pres[(interp_range[2] + 1):length(Tw_time)]) / asc_speed
    }
  }
  
  Ts_star = rep(NA, n_time)
  Ts_star[1] = Tw[1]
  Ts_star_time = Tw_time + delta_t # time associated to Ts* estimates
  
  if ( is.null(MTIME_PARAM) ) {
    #assuming 1m/db, pressure associated to Ts* estimates
    Ts_star_pres = max(Tw_pres) - Ts_star_time*asc_speed
  } else {
    #Ts_star_pres = approx(rev(PRES_PARAM), rev(MTIME_PARAM)*86400, xout=Tw_pres)$y
    Ts_star_pres = approx(rev(MTIME_PARAM)*86400, rev(PRES_PARAM), xout=Ts_star_time)$y
    
    # If necessary, extrapolate the remaining of the axis with 0.1 m/s assumption
    interp_range = range(which(!is.na(Ts_star_pres)))
    if ( interp_range[1] != 1 ) {
      Ts_star_pres[1:(interp_range[1] - 1)] = Ts_star_pres[interp_range[1]] + 
        (Ts_star_time[interp_range[1]] - Ts_star_time[1:(interp_range[1] - 1)]) * asc_speed
    }
    if ( interp_range[2] != length(Ts_star_pres) ) {
      Ts_star_pres[(interp_range[2] + 1):length(Ts_star_pres)] = Ts_star_pres[interp_range[2]] + 
        (Ts_star_time[interp_range[2]] - Ts_star_time[(interp_range[2] + 1):length(Ts_star_pres)]) * asc_speed
    }
  }
  
  for (n in 2:n_time) {
    Ts_star[n] = Ts_star[n-1] + k * (Tw_time[n]-Tw_time[n-1]) * (Tw[n-1]-Ts_star[n-1])
  }
  
  fitted_Ts = approx(Ts_star_pres, Ts_star, xout=PRES_PARAM, rule=2)$y
  
  return(fitted_Ts)
}
