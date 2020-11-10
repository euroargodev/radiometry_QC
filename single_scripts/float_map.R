library(ggplot2)
library(maps)
library(parallel)

world = map_data("world")

corr_list = read.csv("~/Documents/radiometry_QC/float_list.t", sep=" ", header=T)
index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)

#files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
#ident = strsplit(files,"/") #separate the different roots of the files paths
#ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
#dac = ident[,1] #retrieve the DAC of all profiles as a vector
#wod = ident[,2] #retrieve the WMO of all profiles as a vector
#prof_id = ident[,4] #retrieve all profiles  name as a vector
#variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file

corr_list = corr_list[which(corr_list$material!="LIVE"),]

pattern = paste0(corr_list$WMO, "_001.nc")
pattern[126] = paste0(corr_list$WMO[126], "_002.nc") #This one does not have profile 001 ?

n_cores = detectCores()

matches = unlist(mcmapply(grep, pattern, mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(x=index_ifremer$file)))

lat = index_ifremer$latitude[matches]
lon = index_ifremer$longitude[matches]

status = rep(NA, length(lat))
status[which(corr_list$corrected==1)] = "Corrected"
status[which(corr_list$corrected==0)] = "Not corrected"
status[which(is.na(corr_list$corrected))] = "TODO"

lat_plot = lat
lon_plot = lon
buffer = 3 # ?
n = length(lat)
DIST = array(NA, c(n,n))

repeat {
for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		dx = lon_plot[i] - lon_plot[j]
		dy = lat_plot[i] - lat_plot[j]
		d = sqrt(dx^2 + dy^2)
		
		if (d < buffer) {
			my_cos = dx/d
			my_sin = dy/d
			lon_plot[i] = lon_plot[i] + my_cos * (buffer - d)/2
			lon_plot[j] = lon_plot[j] - my_cos * (buffer - d)/2
			lat_plot[i] = lat_plot[i] + my_sin * (buffer - d)/2
			lat_plot[j] = lat_plot[j] - my_sin * (buffer - d)/2
		}
	}
}
for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		dx = lon_plot[i] - lon_plot[j]
		dy = lat_plot[i] - lat_plot[j]
		d = sqrt(dx^2 + dy^2)
		DIST[i,j] = d
	}
}
if (min(DIST, na.rm=T) >= buffer-0.1) {
	break	
}
}


#lat_corr = lat[which(corr_list$corrected==1)]
#lon_corr = lon[which(corr_list$corrected==1)]
#lat_uncorr = lat[which(corr_list$corrected==0)]
#lon_uncorr = lon[which(corr_list$corrected==0)]
#lat_todo = lat[which(is.na(corr_list$corrected))]
#lon_todo = lon[which(is.na(corr_list$corrected))]

g1 = ggplot() +
	geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot, y=lat_plot, color=status)) +
	scale_color_manual(values=c("green", "red", "black"))+
	#geom_point(aes(x=lon_corr, y=lat_corr), color="green") +
	#geom_point(aes(x=lon_uncorr, y=lat_uncorr), color="red") +
	theme_bw()

png(filename="~/Documents/radiometry_QC/float_map.png", width=500, height=300)
g1
dev.off()
