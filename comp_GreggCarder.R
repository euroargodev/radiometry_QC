library(atmos) # Gregg and Carder model
library(ncdf4)
library(stringr)
library(parallel)
library(ggplot2)
library(gridExtra)
library(gtools) # running

n_cores = detectCores()

path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
#path_to_netcdf = "/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)

files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
ident = strsplit(files,"/") #separate the different roots of the files paths
ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
dac = ident[,1] #retrieve the DAC of all profiles as a vector
wod = ident[,2] #retrieve the WMO of all profiles as a vector
prof_id = ident[,4] #retrieve all profiles  name as a vector
variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file
variables = strsplit(variables," ") #separate the different available variables of each profile
lat = index_ifremer$latitude #retrieve the latitude of all profiles as a vector
lon = index_ifremer$longitude #retrieve the longitude of all profiles as a vector
prof_date = index_ifremer$date #retrieve the date of all profiles as a vector


#WMO = "6901473"

### Organelli et al
#WMO = "6901437"
#WMO = "6901510"
#WMO = "6901439"
#WMO = "6901511"
WMO = "6901865"

subset = which(substr(prof_id,3,9)==WMO)

profile_list = paste(path_to_netcdf, files[subset], sep="")

lat_list = lat[subset]
lon_list = lon[subset]

prof_date_list = prof_date[subset]
date_list = as.Date(as.character(prof_date_list), format="%Y%m%d%H%M%S", tz="UTC")
date_list = julian(date_list, origin=as.Date("1950-01-01", tz="UTC"))

month_list = as.numeric(str_sub(prof_date_list,5,6))
day_list = as.numeric(str_sub(prof_date_list,7,8))
hour_list = as.numeric(str_sub(prof_date_list,9,10))
minute_list = as.numeric(str_sub(prof_date_list,11,12))
second_list = as.numeric(str_sub(prof_date_list,13,14))
tu_list = (hour_list + minute_list/60 + second_list/3600)

extrap_Ed <- function(filename) {
    
    filenc = nc_open(filename)
    
    STATION_PARAMETERS = ncvar_get(filenc, "STATION_PARAMETERS")
    param_name = "DOWNWELLING_PAR"

    id_param_arr = which(STATION_PARAMETERS==str_pad(param_name, 64, side="right"), arr.ind=TRUE)
    id_prof = switch (as.character(length(id_param_arr)),
                              "0" = NA,
                              "1" = id_param_arr[2],
                              "2" = id_param_arr[1,2])
    if (is.na(id_prof)) { return(list("DOWN_IRRADIANCE380"=NA, "DOWN_IRRADIANCE412"=NA, "DOWN_IRRADIANCE490"=NA)) }
    n_levels = filenc$dim$N_LEVELS$len
    
    PRES = ncvar_get(filenc, "PRES", start=c(1,id_prof), count=c(n_levels,1))
    
    PARAM_NAMES = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")
    
    all_extrap_Ed = NULL
    
    for (param_name in PARAM_NAMES) {
    
        #IRR = ncvar_get(filenc, param_name)[,id_prof]
        IRR = ncvar_get(filenc, param_name, start=c(1,id_prof), count=c(n_levels,1))
        
        not_na = which(!is.na(IRR) & !is.na(PRES))
        IRR_nna = IRR[not_na]
        PRES_nna = PRES[not_na]
        
        IRR_nna = IRR_nna[order(PRES_nna)]
        PRES_nna = PRES_nna[order(PRES_nna)]
        
        ### find first penetration depth
        
        cut_irr = (1/exp(1)) * max(IRR_nna)
        
        cut_id = max(which(IRR_nna > cut_irr))
        
        if (cut_id < 5) {
            all_extrap_Ed[[param_name]] = NA
            next()
        }
        
        IRR_nna = IRR_nna[1:cut_id]
        PRES_nna = PRES_nna[1:cut_id]
        
        ### second degree polynomial fit
        
        df = data.frame("IRR_nna"=IRR_nna, "PRES_nna"=PRES_nna, "PRES_nna_sq"=PRES_nna^2)
        
        fit_coeff = as.numeric( lm(IRR_nna ~ PRES_nna + PRES_nna_sq, data=df)$coefficients )
        
        all_extrap_Ed[[param_name]] = fit_coeff[1]
        
        #plot(IRR, -PRES)
        #lines(fit_coeff[1] + PRES_nna*fit_coeff[2] + fit_coeff[3]*PRES_nna^2, -PRES_nna, col="red")
    
    }
    
    nc_close(filenc)
    
    return(all_extrap_Ed)
}

valid = which(!is.na(date_list) & !is.na(lon_list) & !is.na(lat_list) & !is.na(tu_list))

Argo_Ed = mcmapply(extrap_Ed, profile_list[valid], mc.cores=n_cores, USE.NAMES=FALSE)

GC_Ed = mcmapply(GreggCarder.f, jday=date_list[valid], rlon=lon_list[valid], rlat=lat_list[valid], 
                 hr=tu_list[valid], MoreArgs=list(lam.sel=c(380,412,490)), mc.cores=n_cores, USE.NAMES=FALSE)[4,]



Argo380 = unlist(Argo_Ed[1,])
Argo412 = unlist(Argo_Ed[2,])
Argo490 = unlist(Argo_Ed[3,])

GC380 = unlist(lapply(GC_Ed, `[[`, 1))
GC412 = unlist(lapply(GC_Ed, `[[`, 2))
GC490 = unlist(lapply(GC_Ed, `[[`, 3))

GC380[which(is.na(GC380))] = 0
GC412[which(is.na(GC412))] = 0
GC490[which(is.na(GC490))] = 0

Argo380_bin = as.vector(running(Argo380, fun=max, width=7, pad=T, by=7))
Argo412_bin = as.vector(running(Argo412, fun=max, width=7, pad=T, by=7))
Argo490_bin = as.vector(running(Argo490, fun=max, width=7, pad=T, by=7))

GC380_bin = as.vector(running(GC380, fun=max, width=7, pad=T, by=7))
GC412_bin = as.vector(running(GC412, fun=max, width=7, pad=T, by=7))
GC490_bin = as.vector(running(GC490, fun=max, width=7, pad=T, by=7))

filter = T
if (filter) {
    Argo380 = as.vector(running(Argo380, fun=max, width=7, pad=T))
    Argo412 = as.vector(running(Argo412, fun=max, width=7, pad=T))
    Argo490 = as.vector(running(Argo490, fun=max, width=7, pad=T))
    
    GC380 = as.vector(running(GC380, fun=max, width=7, pad=T))
    GC412 = as.vector(running(GC412, fun=max, width=7, pad=T))
    GC490 = as.vector(running(GC490, fun=max, width=7, pad=T))
}



PARAM_NAMES = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")

ggdata = data.frame("IRR" = c(Argo380, Argo412, Argo490,
                              GC380, GC412, GC490),
                    "source" = rep(c("Argo","Model"), each=3*length(valid)),
                    "PARAM" = rep(rep(PARAM_NAMES, each=length(valid)), 2),
                    "date" = rep(date_list[valid], 6)
                    )
ggbias = data.frame("IRR" = c(GC380_bin - Argo380_bin, 
                              GC412_bin - Argo412_bin, 
                              GC490_bin - Argo490_bin),
                    "source" = rep("Bias", each=3*length(valid)),
                    "PARAM" = rep(PARAM_NAMES, each=length(valid)),
                    "date" = rep(date_list[valid], 3)
)


g1 = ggplot(ggdata, aes(x=date, y=IRR, color=source, group=source)) +
    geom_line() +
    scale_color_manual(values=c("black", "red")) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(~PARAM, nrow=1)#, scale="free_y")
g2 = ggplot(ggbias, aes(x=date, y=IRR, color=source, group=source)) +
    geom_point(shape=15) +
    scale_color_manual(values=c("blue")) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(~PARAM, nrow=1)#, scale="free_y")

g3 = grid.arrange(g1, g2, nrow=2)

