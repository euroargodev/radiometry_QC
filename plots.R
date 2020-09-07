require(ncdf4)
require(stringr)
require(ggplot2)
require(gridExtra)
require(parallel)
require(tidyr)
require(dplyr)

#require(RColorBrewer)
#require(jsonlite)
#require(networkD3)
#require(tibble)
#require(htmlwidgets)


plot_QC <- function(filename, with_corr=FALSE, do_plot=TRUE, logscale=TRUE) {
 
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
		    nc_close(filenc)
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
		
		if (logscale) {    	
        	g1 = g1 + scale_x_log10()
        	g2 = g2 + scale_x_log10()
        	g3 = g3 + scale_x_log10()
        	g4 = g4 + scale_x_log10()
		} else {
			g1 = g1 + scale_x_continuous(limits=c(-1,1))
			g2 = g2 + scale_x_continuous(limits=c(-1e-4,1e-4))
			g3 = g3 + scale_x_continuous(limits=c(-1e-4,1e-4))
			g4 = g4 + scale_x_continuous(limits=c(-1e-4,1e-4))
		}

		if (i == 1) {
			gg5 = arrangeGrob(g1, g2, g3, g4, nrow=1)
		} else {
			gg5_adjusted = arrangeGrob(g1, g2, g3, g4, nrow=1)
		}
		}
		
	}
    

    if (str_sub(filename,-4,-4) == "D") {
        plot_name = paste("radiometry_QC_", str_sub(filename, -15, -4), "_xing.png", sep="")
    } else {
        plot_name = paste("radiometry_QC_", str_sub(filename, -14, -4), "_xing.png", sep="")
    }

	if (do_plot) {
		if (with_corr){
    		png(plot_name, width = 800, height = 800)
    		plot(grid.arrange(gg5, gg5_adjusted, nrow=2))
    		dev.off()
		} else {
    		png(plot_name, width = 800, height = 400)
			plot(gg5)
    		dev.off()
		}
	}
    
    nc_close(filenc)
    
	return(list("type380"=QC_flags$type380, "type412"=QC_flags$type412, "type490"=QC_flags$type490, "typePAR"=QC_flags$typePAR))

}

plot_corr <- function(filename, pres_zoom=FALSE) {
    filenc = nc_open(filename)
    
    STATION_PARAMETERS = ncvar_get(filenc, "STATION_PARAMETERS")
    param_name = "DOWNWELLING_PAR"
    id_param = grep(str_pad(param_name, 64, side="right"), STATION_PARAMETERS)
    id_prof = arrayInd(id_param, dim(STATION_PARAMETERS))[2]
    if (length(id_param)!=1) {
        nc_close(filenc)
        return(NULL)
    }
		
	PRES = ncvar_get(filenc, "PRES")[,id_prof]
    IRR_380 = ncvar_get(filenc, "DOWN_IRRADIANCE380")[,id_prof]
    IRR_412 = ncvar_get(filenc, "DOWN_IRRADIANCE412")[,id_prof]
    IRR_490 = ncvar_get(filenc, "DOWN_IRRADIANCE490")[,id_prof]
    PAR = ncvar_get(filenc, "DOWNWELLING_PAR")[,id_prof]
    IRR_380_ADJUSTED = ncvar_get(filenc, "DOWN_IRRADIANCE380_ADJUSTED")[,id_prof]
    IRR_412_ADJUSTED = ncvar_get(filenc, "DOWN_IRRADIANCE412_ADJUSTED")[,id_prof]
    IRR_490_ADJUSTED = ncvar_get(filenc, "DOWN_IRRADIANCE490_ADJUSTED")[,id_prof]
    PAR_ADJUSTED = ncvar_get(filenc, "DOWNWELLING_PAR_ADJUSTED")[,id_prof]
	
	nc_close(filenc)

	ggdata = data.frame(PRES, IRR_380, IRR_412, IRR_490, PAR, 
						IRR_380_ADJUSTED, IRR_412_ADJUSTED, IRR_490_ADJUSTED, PAR_ADJUSTED)	
	
	y_limits = NULL
	if (pres_zoom) {
		y_limits = c(250, 0)
	}

	x_limits = NULL

	x_limits[[1]] = c(
		min(min(PAR,na.rm=T), min(PAR_ADJUSTED, na.rm=T)), 
		max(2*quantile(PAR, 1/4, na.rm=T) - min(PAR,na.rm=T), 
			2*quantile(PAR_ADJUSTED, 1/4, na.rm=T) - min(PAR_ADJUSTED, na.rm=T)))
	x_limits[[2]] = c(
		min(min(IRR_380,na.rm=T), min(IRR_380_ADJUSTED, na.rm=T)), 
		max(2*quantile(IRR_380, 1/4, na.rm=T) - min(IRR_380, na.rm=T), 
			2*quantile(IRR_380_ADJUSTED, 1/4, na.rm=T) - min(IRR_380_ADJUSTED, na.rm=T)))
	x_limits[[3]] = c(
		min(min(IRR_412,na.rm=T), min(IRR_412_ADJUSTED, na.rm=T)), 
		max(2*quantile(IRR_412, 1/4, na.rm=T) - min(IRR_412,na.rm=T), 
			2*quantile(IRR_412_ADJUSTED, 1/4, na.rm=T) - min(IRR_412_ADJUSTED, na.rm=T)))
	x_limits[[4]] = c(
		min(min(IRR_490,na.rm=T), min(IRR_490_ADJUSTED, na.rm=T)), 
		max(2*quantile(IRR_490, 1/4, na.rm=T) - min(IRR_490,na.rm=T), 
			2*quantile(IRR_490_ADJUSTED, 1/4, na.rm=T) - min(IRR_490_ADJUSTED, na.rm=T)))
	
    g1 = ggplot(ggdata, aes(x=PAR, y=PRES)) +
        geom_point() +
		geom_path(mapping = aes(x=PAR_ADJUSTED, y=PRES, color="red")) +
       	scale_y_reverse(limits=y_limits) +
       	theme_bw() +
		theme(legend.position="none") 
	
    g2 = ggplot(ggdata, aes(x=IRR_380, y=PRES)) +
        geom_point() +
		geom_path(mapping = aes(x=IRR_380_ADJUSTED, y=PRES, color="red")) +
       	scale_y_reverse(limits=y_limits) +
       	theme_bw() + 
		theme(legend.position="none") 
    
	g3 = ggplot(ggdata, aes(x=IRR_412, y=PRES)) +
        geom_point() +
		geom_path(mapping = aes(x=IRR_412_ADJUSTED, y=PRES, color="red")) +
       	scale_y_reverse(limits=y_limits) +
       	theme_bw() +
		theme(legend.position="none") 
    
	g4 = ggplot(ggdata, aes(x=IRR_490, y=PRES)) +
        geom_point() +
		geom_path(mapping = aes(x=IRR_490_ADJUSTED, y=PRES, color="red")) +
       	scale_y_reverse(limits=y_limits) +
       	theme_bw() +
		theme(legend.position="none") 

    g1_log = g1 + scale_x_log10()
  	g2_log = g2 + scale_x_log10()
   	g3_log = g3 + scale_x_log10()
  	g4_log = g4 + scale_x_log10()

	#g1_lin = g1 + scale_x_continuous(limits=1*c(-1,1))
	#g2_lin = g2 + scale_x_continuous(limits=2*c(-1e-4,1e-4))
	#g3_lin = g3 + scale_x_continuous(limits=3*c(-1e-4,1e-4))
	#g4_lin = g4 + scale_x_continuous(limits=4*c(-1e-4,1e-4))
	g1_lin = g1 + scale_x_continuous(limits=x_limits[[1]])
	g2_lin = g2 + scale_x_continuous(limits=x_limits[[2]])
	g3_lin = g3 + scale_x_continuous(limits=x_limits[[3]])
	g4_lin = g4 + scale_x_continuous(limits=x_limits[[4]])
    
    
	if (str_sub(filename,-4,-4) == "D") {
        plot_name = paste0("radiometry_with_corr_", str_sub(filename, -15, -4), ".png")
    } else {
        plot_name = paste0("radiometry_with_corr_", str_sub(filename, -14, -4), ".png")
    }
    
	png(plot_name, width = 800, height = 800)
	grid.arrange(g1_log, g2_log, g3_log, g4_log, g1_lin, g2_lin, g3_lin, g4_lin, nrow=2)
    dev.off()
	
	return(0)
}

plot_corr_wrapper <- function(WMO, index_ifremer, path_to_netcdf, n_cores=detectCores(), pres_zoom=FALSE) {
    
    files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
    ident = strsplit(files,"/") #separate the different roots of the files paths
    ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
    wod = ident[,2] #retrieve the WMO of all profiles as a vector
    prof_id = ident[,4] #retrieve all profiles  name as a vector
    
    #profile_list = paste(path_to_netcdf, files[which(substr(prof_id,3,9)==WMO)], sep="")
    
    files_corr = files[which(substr(prof_id,3,9)==WMO & substr(prof_id,14,14)==".")]
    for (i in 1:length(files_corr)) {
        path = unlist(strsplit(files_corr[i], "/"))
        str_sub(path[4], 2, 2) = "D"
        files_corr[i] = paste(path[1], path[2], path[3], "RADM/RADM_profiles", path[4], sep="/")
    }
    files_corr = paste(path_to_netcdf, files_corr, sep="")
    
    C = mcmapply(plot_corr, files_corr, mc.cores=n_cores, SIMPLIFY=FALSE, MoreArgs=list(pres_zoom=pres_zoom))
    
    return(C)
}

plot_QC_wrapper <- function(WMO, index_ifremer, path_to_netcdf, n_cores=detectCores(), with_corr=FALSE, do_plot=TRUE, logscale=TRUE) {
    
    files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
    ident = strsplit(files,"/") #separate the different roots of the files paths
    ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
    wod = ident[,2] #retrieve the WMO of all profiles as a vector
    prof_id = ident[,4] #retrieve all profiles  name as a vector
    
    #profile_list = paste(path_to_netcdf, files[which(substr(prof_id,3,9)==WMO)], sep="")
    
    files_corr = files[which(substr(prof_id,3,9)==WMO & substr(prof_id,14,14)==".")]
    for (i in 1:length(files_corr)) {
        path = unlist(strsplit(files_corr[i], "/"))
        files_corr[i] = paste(path[1], path[2], path[3], "RADM/RADM_profiles", path[4], sep="/")
    }
    files_corr = paste(path_to_netcdf, files_corr, sep="")
    
    #M = mcmapply(plot_QC, files_corr, mc.cores=num_cores, SIMPLIFY=FALSE, MoreArgs=list(with_corr=FALSE, do_plot=FALSE))
    
    #Mcorr = mcmapply(plot_QC, files_corr, mc.cores=num_cores, SIMPLIFY=FALSE, MoreArgs=list(with_corr=TRUE, do_plot=FALSE))
    Mcorr = mcmapply(plot_QC, files_corr, mc.cores=num_cores, SIMPLIFY=FALSE, MoreArgs=list(with_corr=with_corr, do_plot=do_plot, logscale=logscale))
    
    return(Mcorr)
    
    
    A = t(array(unlist(M, use.names=FALSE), dim=c(4,length(M))))
    Acorr = t(array(unlist(Mcorr, use.names=FALSE), dim=c(4,length(M))))
    
    param_name = c("IRR380", "IRR412", "IRR490", "PAR")
    sankey_names = paste0("sankey_day_", WMO, "_", param_name, ".html")
    
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
}


#for (file_i in files_corr[260:length(files_corr)]) {
#	apply_plot = plot_QC(file_i, with_corr=TRUE, do_plot=TRUE)
#}
#for (file_i in profile_list[250:270]) {
#	apply_plot = plot_QC(file_i, with_corr=FALSE, do_plot=TRUE)
#}
