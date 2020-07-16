library(atmos) # Gregg and Carder model
library(ncdf4)
library(stringr)
library(parallel)
library(ggplot2)

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


WMO = "6901473"

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

Argo_Ed = mcmapply(extrap_Ed, profile_list, mc.cores=n_cores, USE.NAMES=FALSE)

GC_Ed = mcmapply(GreggCarder.f, jday=date_list, rlon=lon_list, rlat=lat_list, hr=tu_list, MoreArgs=list(lam.sel=c(380,412,490)), 
                 mc.cores=n_cores, USE.NAMES=FALSE)[4,]

PARAM_NAMES = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", "DOWN_IRRADIANCE490")

ggdata = data.frame("IRR" = c(unlist(Argo_Ed[1,]), unlist(Argo_Ed[2,]), unlist(Argo_Ed[3,]),
                              unlist(lapply(GC_Ed, `[[`, 1)), unlist(lapply(GC_Ed, `[[`, 2)), unlist(lapply(GC_Ed, `[[`, 3)),
                              (unlist(lapply(GC_Ed, `[[`, 1)) - unlist(Argo_Ed[1,])), 
                              (unlist(lapply(GC_Ed, `[[`, 2)) - unlist(Argo_Ed[2,])), 
                              (unlist(lapply(GC_Ed, `[[`, 3)) - unlist(Argo_Ed[3,]))),
                    "source" = rep(c("Argo","Model","Bias"), each=3*length(profile_list)),
                    "PARAM" = rep(rep(PARAM_NAMES, each=length(profile_list)), 3),
                    "date" = rep(date_list, 9)
                    )

g1 = ggplot(na.omit(ggdata), aes(x=date, y=IRR, color=source, group=PARAM)) +
    geom_point() + #data=function(x){x[!x$is_greylisted, ]}) +
    #geom_point(data=function(x){x[x$is_greylisted, ]}, color="red") +
    #scale_color_viridis() +
    #scale_y_continuous(
    #    name = "First Axis",
    #    sec.axis = sec_axis(~./2.5, name="Second Axis")
    #) +
    facet_wrap(~PARAM, ncol=1)#, scale="free_y")
