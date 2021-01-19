require(ncdf4)
require(stringr)

df = read.table("~/Documents/radiometry_QC/floats_DM_radiometry.txt", col.names=c("WMO"))

PARAM_NAMES = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", 
				"DOWN_IRRADIANCE490", "DOWNWELLING_PAR")

for (PARAM in PARAM_NAMES) {
	df[paste0(PARAM, "_A")] = rep(NA, length(df$WMO))
	df[paste0(PARAM, "_B")] = rep(NA, length(df$WMO))
	df[paste0(PARAM, "_C")] = rep(NA, length(df$WMO))
	df[paste0(PARAM, "_Q")] = rep(NA, length(df$WMO))
}

for (i in 1:length(df$WMO)) {
	for (PARAM in PARAM_NAMES) {
		filename = paste0("/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/coriolis/", 
					  	  df$WMO[i], "/profiles/BD", df$WMO[i], "_001.nc")
		fnc = nc_open(filename)
		coef = ncvar_get(fnc, "SCIENTIFIC_CALIB_COEFFICIENT")	
		param = ncvar_get(fnc, "PARAMETER")
		id = which(param == str_pad(PARAM, 64, side="right"))
		calib = unlist(str_split(coef[id], ","))
		for (letter in c("A", "B", "C", "Q")) {
			if (any(grepl(letter, calib))) {
				value = unlist(str_split(calib[grepl(letter, calib)], "="))[2]
				df[i, paste0(PARAM, "_", letter)] = as.numeric(value)
			} else {
				df[i, paste0(PARAM, "_", letter)] = 0
			}
		}
	}
}
