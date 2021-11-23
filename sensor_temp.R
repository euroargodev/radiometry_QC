sensor_temp <- function(TEMP, PRES_TEMP, PRES_PARAM, material="PEEK", MTIME_PARAM) {
  
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
  
  #asc_speed = 0.1 # m/s

  Tw = rev(TEMP[!is.na(TEMP) & !is.na(PRES_TEMP)])

  #pressure associated to Tw measurements
  Tw_pres = rev(PRES_TEMP[!is.na(TEMP) & !is.na(PRES_TEMP)])
  
  ##Time associated to Tw measurements
  #Tw_time = (max(Tw_pres) - Tw_pres) / asc_speed

  n_time = length(Tw)
  
  if (n_time < 2) {
    return(rep(NA, length(PRES_PARAM)))
  }
  
  flag_PARAM = !is.na(PRES_PARAM)
  
  Tw_time = (rev(MTIME_PARAM[flag_PARAM])-min(MTIME_PARAM[flag_PARAM]))*86400
  
  Tw_PARAM = approx(Tw_pres, Tw, xout=rev(PRES_PARAM[flag_PARAM]),rule=2)$y
  
  Tw_star_PARAM = approx(Tw_time, Tw_PARAM, xout=Tw_time - delta_t,rule=2)$y
  
  n_time = length(Tw_PARAM)
  
  Ts = rep(NA, n_time)
  Ts[1] = Tw_star_PARAM[1]
  #Ts_time = Tw_time + delta_t # time associated to Ts estimates
  #Ts_pres = max(Tw_pres) - Ts_time*asc_speed # pressure associated to Ts estimates

  for (n in 2:n_time) {
    Ts[n] = Ts[n-1] + k * (Tw_time[n]-Tw_time[n-1]) * (Tw_star_PARAM[n-1]-Ts[n-1])
  }
  
  #fitted_Ts = approx(Ts_pres, Ts, xout=PRES_PARAM)$y
  #fitted_Ts = rev(Ts)
  fitted_Ts = PRES_PARAM*NA
  fitted_Ts[flag_PARAM] = rev(Ts)
  
  return(fitted_Ts)
}
