library(ncdf4)
library(stringr)
library(ggplot2)
library(gridExtra)
library(latex2exp)

source("~/Documents/radiometry/get_matches.R")
source("~/Documents/radiometry/sensor_temp.R")

WMO = "6901474"
profile = "006"
#param = "DOWN_IRRADIANCE380"
#param = "DOWNWELLING_PAR"
PARAM_NAMES = c("DOWNWELLING_PAR", "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")
irr_labels = c("PAR", "E_d(380)", "E_d(412)", "E_d(490)")
names(irr_labels) = PARAM_NAMES
units = c("$\\mu mol\\cdot m^{-2}\\cdot s^{-1}$", "$W\\cdot m^{-2}\\cdot nm^{-1}$", "$W\\cdot m^{-2}\\cdot nm^{-1}$", "$W\\cdot m^{-2}\\cdot nm^{-1}$")
names(units) = PARAM_NAMES
limits = list()
limits[[PARAM_NAMES[1]]] = c(-1e-1, 3e-1)
limits[[PARAM_NAMES[2]]] = c(-1e-4, 3e-4)
limits[[PARAM_NAMES[3]]] = c(-1e-4, 3e-4)
limits[[PARAM_NAMES[4]]] = c(-1e-4, 3e-4)


file_A = paste0("/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/coriolis/", WMO, "/profiles/BD", WMO, "_", profile, ".nc")
file_D = paste0("~/Documents/radiometry_QC/", WMO, "/RADM_profiles/BD", WMO, "_", profile, ".nc")
file_core = paste0("/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/coriolis/", WMO, "/profiles/D", WMO, "_", profile, ".nc")

path_to_netcdf=""

get_param <- function(filename, PARAM_NAME, adj = T) {
    fnc = nc_open(filename)
    
    parameters = ncvar_get(fnc, "STATION_PARAMETERS")
    id_param_arr = which(parameters == str_pad(PARAM_NAME, 64, side="right"), arr.ind=TRUE)
    if (length(id_param_arr)==2) id_prof=id_param_arr[2] else id_prof=id_param_arr[1,2]
    if (length(id_param_arr)==2) id_param=id_param_arr[1] else id_param=id_param_arr[1,1]
    
    n_level = fnc$dim$N_LEVEL$len
    
    PARAM_NAME_ADJUSTED = paste0(PARAM_NAME, "_ADJUSTED")
    
    if (adj) {
        param_return = ncvar_get(fnc, PARAM_NAME_ADJUSTED, start=c(1,id_prof), count=c(n_level,1))	
    } else {
        param_return = ncvar_get(fnc, PARAM_NAME, start=c(1,id_prof), count=c(n_level,1))
    }
    nc_close(fnc)
    
    return(param_return)
}

for (param in PARAM_NAMES) { 

match_A = get_Ts_match(path_to_netcdf, file_A, param)
Irr_D = get_param(file_D, param)

nna = which(!is.na(match_A$PARAM) & !is.na(match_A$Ts) & !is.na(match_A$PRES))

df = data.frame("PARAM_A"=match_A$PARAM[nna], "PARAM_D"=Irr_D[nna], "PRES"=match_A$PRES[nna], "Ts"=match_A$Ts[nna])
df$ERROR = pmax(2.5e-5, df$PARAM_D*0.02)

ab = lm(df$PARAM_D-df$PARAM_A ~ df$Ts)

g1 = ggplot(df, aes(x=Ts, y=PRES)) +
	scale_y_reverse(name="Pressure (dbar)") +
	scale_x_continuous(
		name = "Sensor Temperature (°C)",
		sec.axis = sec_axis(trans=~.*ab$coefficients[2] + ab$coefficients[1], name=TeX(paste0("Applied offset (", units[param], ")")), 
		#breaks = waiver()
		#breaks = c(-4e-5, -1.5e-5, 1e-5)
		)
	) +
	#geom_ribbon(aes(xmin=Ts-2.5, xmax=Ts+2.5), fill="grey70") +
	geom_path() +
	theme_bw() +
	theme(plot.margin=margin(t=5.5, b=5.5, l=5.5, r=15, unit="pt"))

g2 = ggplot(df, aes(x=PARAM_A, y=PRES)) +
	scale_y_reverse(name="Pressure (dbar)") +
	scale_x_continuous(
		limits=limits[[param]],
		name = TeX(paste0(irr_labels[param]," (", units[param], ")")),
		sec.axis = sec_axis(trans=~.*1, name=TeX(paste0("Applied offset (", units[param], ")"))) 
) +
	geom_vline(xintercept=0, linetype="dashed") +
	scale_color_manual(values = c("Not corrected"="black", "Corrected"="red")) +
	geom_ribbon(aes(xmin=PARAM_D-ERROR, xmax=PARAM_D+ERROR), fill="grey60", alpha=0.5) +
	geom_path(aes(color="Not corrected")) +
	geom_path(aes(x=PARAM_D, color="Corrected")) +
	theme_bw() +
	theme(
		axis.title.x.top = element_text(color = "white"),
		axis.text.x.top = element_text(color = "white"),
		axis.ticks.x.top = element_line(size=0),
		legend.title = element_blank(),
		legend.position = c(0.01, 0.99),
		legend.justification = c("left", "top"),
		legend.background = element_rect(fill="grey90")
	)

min_ribbon = min(df$PARAM_A[which(df$PARAM_A>0)], df$PARAM_D[which(df$PARAM_D>0)])

g3 = ggplot(df, aes(x=PARAM_A, y=PRES)) +
	scale_y_reverse(name="Pressure (dbar)") +
	scale_x_log10(
		expand=expansion(mult=c(0,0.04)),
		name = TeX(paste0(irr_labels[param]," (", units[param], ")")),
		sec.axis = sec_axis(trans=~.*1, name=TeX(paste0("Applied offset (", units[param], ")"))) 
) + 
	#geom_vline(xintercept=0, linetype="dashed") +
    geom_ribbon(aes(xmin=pmax(min_ribbon,PARAM_D-ERROR), xmax=PARAM_D+ERROR), fill="grey60", alpha=0.5) +
	geom_point(aes(color="Not corrected"), size=1) +
	geom_point(aes(x=PARAM_D, color="Corrected"), size=1) +
	scale_color_manual(values = c("Not corrected"="black", "Corrected"="red")) +
	theme_bw() +
	theme(
		axis.title.x.top = element_text(color = "white"),
		axis.text.x.top = element_text(color = "white"),
		axis.ticks.x.top = element_line(size=0),
		legend.title = element_blank(),
		legend.position = c(0.01, 0.99),
		legend.justification = c("left", "top"),
		legend.background = element_rect(fill="grey90")
	) 

#gt = grid.arrange(g3, g2, g1, nrow=1)
#gt

png(filename=paste0("~/Documents/radiometry_QC/plot_diff/plot_correction_", WMO, "_", profile, "_", param, ".png"), width=800, height=300)
grid.arrange(g3, g2, g1, nrow=1)
dev.off()

}
