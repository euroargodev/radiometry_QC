library(ggplot2)
library(maps)
library(parallel)
library(stringr)
library(wesanderson)
library(RColorBrewer)

world = map_data("world")

corr_list = read.csv("~/Documents/radiometry_QC/float_list_temp.t", sep=" ", header=T)
index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)

#files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
#ident = strsplit(files,"/") #separate the different roots of the files paths
#ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
#dac = ident[,1] #retrieve the DAC of all profiles as a vector
#wod = ident[,2] #retrieve the WMO of all profiles as a vector
#prof_id = ident[,4] #retrieve all profiles  name as a vector
#variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file

corr_list = corr_list[which(corr_list$material!="LIVE"),]

pattern = paste0(corr_list$WMO, "_[0-9][05]1.nc")
pattern[126] = paste0(corr_list$WMO[126], "_[0-9][05]2.nc") #This one does not have profile 001 ?

n_cores = detectCores()

matches = unlist(mcmapply(grep, pattern, mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(x=index_ifremer$file)))

lat = index_ifremer$latitude[matches]
lon = index_ifremer$longitude[matches]

match_WMO = str_sub(index_ifremer$file[matches], 10, 16)

status = rep(NA, length(lat))
status[which(corr_list$corrected[match(match_WMO, corr_list$WMO)]==1)] = "Corrected"
status[which(corr_list$corrected[match(match_WMO, corr_list$WMO)]==0)] = "Not corrected"
#status[which(is.na(corr_list$corrected))] = "TODO"

lat_range = c(30, 47)
lon_range = c(-6, 43)
medblack = ( lat > lat_range[1] & lat < lat_range[2] & lon > lon_range[1] & lon < lon_range[2])

lat_plot = lat[which(!medblack)]
lon_plot = lon[which(!medblack)]
status_plot = status[which(!medblack)]
buffer = 2 # ?
n = length(lat_plot)
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

lat_plot_mb = lat[which(medblack)]
lon_plot_mb = lon[which(medblack)]
status_plot_mb = status[which(medblack)]
buffer_mb = 0.8 # ?
n_mb = length(lat_plot_mb)
DIST_mb = array(NA, c(n,n))
#map_mb = map("world", fill=TRUE, xlim = lon_range, ylim = lat_range)
subworld_mb = which(world$long > lon_range[1] & world$long < lon_range[2] &
					world$lat > lat_range[1] & world$lat < lat_range[2])
map_mb = world[subworld_mb,]

repeat {
for (i in 1:(n_mb-1)) {
	for (j in (i+1):n_mb) {
		dx = lon_plot_mb[i] - lon_plot_mb[j]
		dy = lat_plot_mb[i] - lat_plot_mb[j]
		d = sqrt(dx^2 + dy^2)
		if (d < buffer_mb) {
			my_cos = dx/d
			my_sin = dy/d
			lon_plot_mb[i] = lon_plot_mb[i] + my_cos * (buffer_mb - d)/2
			lon_plot_mb[j] = lon_plot_mb[j] - my_cos * (buffer_mb - d)/2
			lat_plot_mb[i] = lat_plot_mb[i] + my_sin * (buffer_mb - d)/2
			lat_plot_mb[j] = lat_plot_mb[j] - my_sin * (buffer_mb - d)/2
		}
	}
}
for (i in 1:(n_mb-1)) {
	for (j in (i+1):n_mb) {
		dx = lon_plot_mb[i] - lon_plot_mb[j]
		dy = lat_plot_mb[i] - lat_plot_mb[j]
		d = sqrt(dx^2 + dy^2)
		DIST_mb[i,j] = d
	}
}
if (min(DIST_mb, na.rm=T) >= buffer_mb-0.1) {
	break	
}
}


#my_palette = c(wes_palette("Rushmore1")[3], wes_palette("Rushmore1")[5])
#my_palette = c(wes_palette("Cavalcanti1")[2], wes_palette("Cavalcanti1")[5])
my_palette = brewer.pal(n=3, name="Dark2")[1:2]


g1 = ggplot() +
	geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot, y=lat_plot, color=status_plot), size=0.7) +
	scale_color_manual(values=my_palette)+
	annotate("rect", xmin=lon_range[1], xmax=lon_range[2], ymin=lat_range[1], ymax = lat_range[2], color="black", fill="black", alpha=0.5) +
	theme_bw() +
	labs(color = "Correction status", x="Longitude", y="Latitude") 
g2 = ggplot() +
	geom_polygon(data=map_mb, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot_mb, y=lat_plot_mb, color=status_plot_mb), size=0.7) +
	scale_color_manual(values=my_palette)+
	theme_bw() +
	theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

png(filename="~/Documents/radiometry_QC/profile_map.png", width=700, height=400)
plot(g1)
dev.off()
png(filename="~/Documents/radiometry_QC/profile_map_medblack.png", width=200, height=150)
plot(g2)
dev.off()
