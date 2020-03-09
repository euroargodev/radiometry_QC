library(ncdf4)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

source("RT_QC_radiometry_function_oao_2.R")

path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

index_ifremer = read.table("argo_bio-profile_index.txt", sep=",", header = T)

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


profile_actual="6901524_150."

iii = which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index

filename = paste(path_to_netcdf, files[iii], sep="")


filenc = nc_open(filename)

STATION_PARAMETERS = ncvar_get(filenc, "STATION_PARAMETERS")
param_name = "DOWNWELLING_PAR"
id_prof = which(STATION_PARAMETERS==str_pad(param_name, 64, side="right"), arr.ind=TRUE)[2]

PRES = ncvar_get(filenc, "PRES")[,id_prof]
IRR_380 = ncvar_get(filenc, "DOWN_IRRADIANCE380")[,id_prof]
IRR_412 = ncvar_get(filenc, "DOWN_IRRADIANCE412")[,id_prof]
IRR_490 = ncvar_get(filenc, "DOWN_IRRADIANCE490")[,id_prof]
PAR = ncvar_get(filenc, "DOWNWELLING_PAR")[,id_prof]

QC_flags = RT_QC_radiometry(PRES, IRR_380, IRR_412, IRR_490, PAR)

dataf = data.frame(QC_flags)
dataf$PRES = PRES
dataf$IRR_380 = IRR_380
dataf$IRR_412 = IRR_412
dataf$IRR_490 = IRR_490
dataf$PAR = PAR

dataf$FLAG_380[which(dataf$FLAG_380 == "NA")] = NA
dataf$FLAG_412[which(dataf$FLAG_412 == "NA")] = NA
dataf$FLAG_490[which(dataf$FLAG_490 == "NA")] = NA
dataf$FLAG_PAR[which(dataf$FLAG_PAR == "NA")] = NA

dataf$FLAG_PAR = factor(dataf$FLAG_PAR , levels=c("1", "2", "3", "4", NA))

col_1 = "#869C66"
col_2 = "#FFCA0B"
col_3 = "#FF9502"
col_4 = "#D43501"

type_col = c(col_1, col_2, col_3)

g1 = ggplot(dataf, aes(y=PRES, x=PAR, color=FLAG_PAR)) +
    geom_point() +
    scale_y_reverse() +
    scale_x_log10() +
    scale_color_manual(values=c(col_1, col_2, col_3, col_4)) +
    geom_label(
        label=paste("Profile type :", dataf$typePAR), 
        x=log10(max(dataf$PAR, na.rm=T)),
        y=-max(dataf$PRES, na.rm=T),
        vjust="bottom",
        hjust="right",
        label.padding = unit(0.55, "lines"), # Rectangle size around label
        #label.size = 0.35,
        color = "black",
        fill=type_col[as.numeric(dataf$typePAR)]
    ) +
    theme_bw() +
    theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))

g2 = ggplot(dataf, aes(y=PRES, x=IRR_380, color=FLAG_380)) +
    geom_point() +
    scale_y_reverse() +
    scale_x_log10() +
    scale_color_manual(values=c(col_1, col_2, col_3, col_4)) +
    geom_label(
        label=paste("Profile type :", dataf$type380), 
        x=log10(max(dataf$IRR_380, na.rm=T)),
        y=-max(dataf$PRES, na.rm=T),
        vjust="bottom",
        hjust="right",
        label.padding = unit(0.55, "lines"), # Rectangle size around label
        #label.size = 0.35,
        color = "black",
        fill=type_col[as.numeric(dataf$type380)]
    ) +
    theme_bw() +
    theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))

g3 = ggplot(dataf, aes(y=PRES, x=IRR_412, color=FLAG_412)) +
    geom_point() +
    scale_y_reverse() +
    scale_x_log10() +
    scale_color_manual(values=c(col_1, col_2, col_3, col_4)) +
    geom_label(
        label=paste("Profile type :", dataf$type412), 
        x=log10(max(dataf$IRR_412, na.rm=T)),
        y=-max(dataf$PRES, na.rm=T),
        vjust="bottom",
        hjust="right",
        label.padding = unit(0.55, "lines"), # Rectangle size around label
        #label.size = 0.35,
        color = "black",
        fill=type_col[as.numeric(dataf$type412)]
    ) +
    theme_bw() +
    theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))

g4 = ggplot(dataf, aes(y=PRES, x=IRR_490, color=FLAG_490)) +
    geom_point() +
    scale_y_reverse() +
    scale_x_log10() +
    scale_color_manual(values=c(col_1, col_2, col_3, col_4)) +
    geom_label(
        label=paste("Profile type :", dataf$type490), 
        x=log10(max(dataf$IRR_490, na.rm=T)),
        y=-max(dataf$PRES, na.rm=T),
        vjust="bottom",
        hjust="right",
        label.padding = unit(0.55, "lines"), # Rectangle size around label
        #label.size = 0.35,
        color = "black",
        fill=type_col[as.numeric(dataf$type490)]
    ) +
    theme_bw() +
    theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))


grid.arrange(g1, g2, g3, g4, nrow=1)



