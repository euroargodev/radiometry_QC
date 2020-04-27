library(ncdf4)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(parallel)
library(jsonlite)
library(networkD3)
library(dplyr)
library(tibble)
library(tidyr)
library(htmlwidgets)

source("~/Documents/radiometry/RT_QC_radiometry_function_oao_2.R")

#path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
path_to_netcdf = "/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

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


WMO = "6901525"

profile_list = paste(path_to_netcdf, files[which(substr(prof_id,3,9)==WMO)], sep="")

files_corr = files[which(substr(prof_id,3,9)==WMO & substr(prof_id,14,14)==".")]
for (i in 1:length(files_corr)) {
	path = unlist(strsplit(files_corr[i], "/"))
	files_corr[i] = paste(path[1], path[2], path[3], "radiometry_xing", path[4], sep="/")
}
files_corr = paste(path_to_netcdf, files_corr, sep="")


plot_QC <- function(filename, with_corr=FALSE, do_plot=TRUE) {
 
    filenc = nc_open(filename)
    
    STATION_PARAMETERS = ncvar_get(filenc, "STATION_PARAMETERS")
    param_name = "DOWNWELLING_PAR"
    id_param = grep(str_pad(param_name, 64, side="right"), STATION_PARAMETERS)
    id_prof = arrayInd(id_param, dim(STATION_PARAMETERS))[2]
    if (length(id_param)!=1) {
        nc_close(filenc)
        return(NULL)
    }

	IRR380_NAME = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE380_ADJUSTED")
	IRR412_NAME = c("DOWN_IRRADIANCE412", "DOWN_IRRADIANCE412_ADJUSTED")
	IRR490_NAME = c("DOWN_IRRADIANCE490", "DOWN_IRRADIANCE490_ADJUSTED")
	PAR_NAME = c("DOWNWELLING_PAR", "DOWNWELLING_PAR_ADJUSTED")
    
	if (with_corr) {
		max_i = 2
	} else {
		max_i = 1
	}	

	for (i in 1:max_i) {
    
		PRES = ncvar_get(filenc, "PRES")[,id_prof]
    	IRR_380 = ncvar_get(filenc, IRR380_NAME[i])[,id_prof]
    	IRR_412 = ncvar_get(filenc, IRR412_NAME[i])[,id_prof]
    	IRR_490 = ncvar_get(filenc, IRR490_NAME[i])[,id_prof]
    	PAR = ncvar_get(filenc, PAR_NAME[i])[,id_prof]

		if (length(which(!is.na(IRR_380))) < 5) {
			return(list("type380"="0", "type412"="0", "type490"="0", "typePAR"="0"))
		}   
 
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
    
    	dataf$FLAG_380 = factor(dataf$FLAG_380 , levels=c("1", "2", "3", "4", NA))
    	dataf$FLAG_412 = factor(dataf$FLAG_412 , levels=c("1", "2", "3", "4", NA))
    	dataf$FLAG_490 = factor(dataf$FLAG_490 , levels=c("1", "2", "3", "4", NA))
    	dataf$FLAG_PAR = factor(dataf$FLAG_PAR , levels=c("1", "2", "3", "4", NA))
   		
		if (do_plot) { 
    	col_1 = "#869C66"
    	col_2 = "#FFCA0B"
    	col_3 = "#FF9502"
    	col_4 = "#D43501"
    	
    	type_col = c(col_1, col_2, col_3)
    
    	g1 = ggplot(dataf, aes(y=PRES, x=PAR, color=FLAG_PAR)) +
    	    geom_point() +
        	scale_y_reverse() +
        	scale_x_log10() +
        	scale_color_manual(values=c(col_1, col_2, col_3, col_4), breaks=c("1", "2", "3", "4"), drop=FALSE) +
        	geom_label(
            	label=paste("Profile type :", dataf$typePAR), 
            	x=log10(max(dataf$PAR, na.rm=T)),
            	y=-max(dataf$PRES, na.rm=T),
            	vjust="bottom",
            	hjust="right",
            	label.padding = unit(0.55, "lines"), # Rectangle size around label
            	#label.size = 0.35,
            	color = "black",
            	fill = type_col[as.numeric(as.character(dataf$typePAR[1]))]
        	) +
        	theme_bw() +
        	theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))
    
   		g2 = ggplot(dataf, aes(y=PRES, x=IRR_380, color=FLAG_380)) +
	        geom_point() +
	        scale_y_reverse() +
	        scale_x_log10() +
	        scale_color_manual(values=c(col_1, col_2, col_3, col_4), breaks=c("1", "2", "3", "4"), drop=FALSE) +
	        geom_label(
	            label=paste("Profile type :", dataf$type380), 
	            x=log10(max(dataf$IRR_380, na.rm=T)),
	            y=-max(dataf$PRES, na.rm=T),
	            vjust="bottom",
	            hjust="right",
	            label.padding = unit(0.55, "lines"), # Rectangle size around label
	            #label.size = 0.35,
	            color = "black",
	            fill=type_col[as.numeric(as.character(dataf$type380[1]))]
	        ) +
	        theme_bw() +
	        theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))
	    
	    g3 = ggplot(dataf, aes(y=PRES, x=IRR_412, color=FLAG_412)) +
	        geom_point() +
	        scale_y_reverse() +
	        scale_x_log10() +
	        scale_color_manual(values=c(col_1, col_2, col_3, col_4), breaks=c("1", "2", "3", "4"), drop=FALSE) +
	        geom_label(
	            label=paste("Profile type :", dataf$type412), 
	            x=log10(max(dataf$IRR_412, na.rm=T)),
	            y=-max(dataf$PRES, na.rm=T),
	            vjust="bottom",
	            hjust="right",
	            label.padding = unit(0.55, "lines"), # Rectangle size around label
	            #label.size = 0.35,
	            color = "black",
	            fill=type_col[as.numeric(as.character(dataf$type412[1]))]
	        ) +
	        theme_bw() +
	        theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))
	    
	    g4 = ggplot(dataf, aes(y=PRES, x=IRR_490, color=FLAG_490)) +
	        geom_point() +
	        scale_y_reverse() +
	        scale_x_log10() +
	        scale_color_manual(values=c(col_1, col_2, col_3, col_4), breaks=c("1", "2", "3", "4"), drop=FALSE) +
	        geom_label(
	   	        label=paste("Profile type :", dataf$type490), 
	            x=log10(max(dataf$IRR_490, na.rm=T)),
	            y=-max(dataf$PRES, na.rm=T),
	            vjust="bottom",
	            hjust="right",
	            label.padding = unit(0.55, "lines"), # Rectangle size around label
	            #label.size = 0.35,
	            color = "black",
	            fill=type_col[as.numeric(as.character(dataf$type490[1]))]
	        ) +
	        theme_bw() +
	        theme(legend.justification=c("right", "bottom"), legend.position=c(0.95,0.15))
    	
		if (i == 1) {
			gg5 = grid.arrange(g1, g2, g3, g4, nrow=1)
		} else {
			gg5_adjusted = grid.arrange(g1, g2, g3, g4, nrow=1)
		}
		}
		
	}
    

    if (str_sub(filename,-4,-4) == "D") {
        plot_name = paste("radiometry_QC_", str_sub(filename, -15, -4), ".png", sep="")
    } else {
        plot_name = paste("radiometry_QC_", str_sub(filename, -14, -4), ".png", sep="")
    }

	if (do_plot) {
		if (with_corr){
    		png(plot_name, width = 800, height = 800)
    		grid.arrange(gg5, gg5_adjusted, nrow=2)
    		dev.off()
		} else {
    		png(plot_name, width = 800, height = 400)
			grid.arrange(gg5)
    		dev.off()
		}
	}
    
    nc_close(filenc)
    
	return(list("type380"=QC_flags$type380, "type412"=QC_flags$type412, "type490"=QC_flags$type490, "typePAR"=QC_flags$typePAR))

}


num_cores = detectCores()
M = mcmapply(plot_QC, files_corr, mc.cores=num_cores, SIMPLIFY=FALSE, MoreArgs=list(with_corr=FALSE, do_plot=FALSE))

Mcorr = mcmapply(plot_QC, files_corr, mc.cores=num_cores, SIMPLIFY=FALSE, MoreArgs=list(with_corr=TRUE, do_plot=FALSE))

A = t(array(unlist(M, use.names=FALSE), dim=c(4,length(M))))
Acorr = t(array(unlist(Mcorr, use.names=FALSE), dim=c(4,length(M))))

param_name = c("IRR380", "IRR412", "IRR490", "PAR")
sankey_names = paste0("sankey_", WMO, "_", param_name, ".html")

for (n in 1:4) {
	type = A[,n]
	type_corr = Acorr[,n]

	change = matrix(0, 4, 4)
	rownames(change) = c("1", "2", "3", "4")
	colnames(change) = c("corr_1", "corr_2", "corr_3", "corr_4")

	for (i in 1:4) {
		for (j in 1:4) {
			change[i,j] = length(which(type==as.character(i) & type_corr==as.character(j)))
		}
	} 

	links = change %>%
		as.data.frame() %>%
		rownames_to_column(var="source") %>%
		gather(key="target", value="value", -1) %>%
		filter(value != 0)

	nodes = data.frame(
		name = c(as.character(links$source), as.character(links$target)) %>%
		unique()
		)

	links$IDsource = match(links$source, nodes$name) -1
	links$IDtarget = match(links$target, nodes$name) -1

	my_color = 'd3.scaleOrdinal() .domain(["1", "2", "3", "4", "corr_1", "corr_2", "corr_3", "corr_4"]) .range(["#869C66", "#FFCA0B", "#FF9502", "#D43501", "#869C66", "#FFCA0B", "#FF9502", "#D43501"])'

	p1 = sankeyNetwork(Links = links, Nodes = nodes,
					Source = "IDsource", Target = "IDtarget",
					Value = "value", NodeID = "name",
					colourScale = my_color)

	saveWidget(p1, file=sankey_names[n], selfcontained=FALSE)
}


#write_json(M, path="profile_types.json", auto_unbox=T, pretty=T, null="list")
