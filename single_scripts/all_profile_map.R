library(ggplot2)
library(maps)
library(parallel)
library(stringr)
library(wesanderson)
library(RColorBrewer)

world = map_data("world")

#corr_list = read.csv("~/Documents/radiometry_QC/float_list_temp.t", sep=" ", header=T)
#index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)
source("~/Documents/radiometry/single_scripts/survey.R")

#corr_list = corr_list[which(corr_list$material!="LIVE"),]
corr_list = corr_list[which(corr_list$material!="LIVE" & corr_list$num_night>0 & corr_list$num_drift>=0.8*corr_list$num_prof),]
#corrected_with_alternative = c("6901580", "6901866", "6901650", "6902968", "6901647", "6902547") # these floats have good QC data but were corrected with B and/or day methods
corrected_with_alternative = c("6901580", "6901866", "6901650", "6902968") # these floats have good QC data but were corrected with B and/or day methods (removed the 2 floats that could be corrected with just main methods)
corr_list$corrected_with_alternative = rep(FALSE, length(corr_list$WMO))
for (the_WMO in corrected_with_alternative) {
	corr_list$corrected_with_alternative[which(corr_list$WMO == the_WMO)] = TRUE
}


pattern = paste0(corr_list$WMO)
#pattern[126] = paste0(corr_list$WMO[126], "_[0-9][05]2.nc") #This one does not have profile 001 ?

n_cores = detectCores()

matches = unlist(mcmapply(grep, pattern, mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(x=index_ifremer$file)))

lat = index_ifremer$latitude[matches]
lon = index_ifremer$longitude[matches]

match_WMO = str_sub(index_ifremer$file[matches], 10, 16)

status = rep(NA, length(lat))
status[which(corr_list$corrected[match(match_WMO, corr_list$WMO)]==1)] = "Corrected"
status[which(corr_list$corrected[match(match_WMO, corr_list$WMO)]==0)] = "Not corrected"
status[which(corr_list$corrected_with_alternative[match(match_WMO, corr_list$WMO)])] = "Corrected with alternative method(s)"
#status[which(is.na(corr_list$corrected))] = "TODO"

lat_range = c(30, 47)
lon_range = c(-6, 43)
medblack = ( lat > lat_range[1] & lat < lat_range[2] & lon > lon_range[1] & lon < lon_range[2])

lat_plot = lat[which(!medblack)]
lon_plot = lon[which(!medblack)]
status_plot = status[which(!medblack)]


lat_plot_mb = lat[which(medblack)]
lon_plot_mb = lon[which(medblack)]
status_plot_mb = status[which(medblack)]

subworld_mb = which(world$long > lon_range[1] & world$long < lon_range[2] &
					world$lat > lat_range[1] & world$lat < lat_range[2])
map_mb = world[subworld_mb,]
map_world = world[which(world$long >= -180 & world$long <= 180),]


#my_palette = c(wes_palette("Rushmore1")[3], wes_palette("Rushmore1")[5])
#my_palette = c(wes_palette("Cavalcanti1")[2], wes_palette("Cavalcanti1")[5])
my_palette = brewer.pal(n=6, name="Dark2")[c(1,2,6,3)]


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

g3 = ggplot() +
	geom_polygon(data=map_world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot[status_plot=="Corrected"], y=lat_plot[status_plot=="Corrected"]), color=my_palette[1], size=0.7) +
	annotate("rect", xmin=lon_range[1], xmax=lon_range[2], ymin=lat_range[1], ymax = lat_range[2], color="black", fill="black", alpha=0.5) +
	scale_x_continuous(breaks=seq(-180,180,30)) +
	scale_y_continuous(breaks=seq(-90,90,30)) +
	theme_bw() +
	labs(color = "Correction status", x="Longitude", y="Latitude") 
g3.5 = ggplot() +
	geom_polygon(data=map_world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot[status_plot=="Corrected"], y=lat_plot[status_plot=="Corrected"]), color=my_palette[1], size=0.7) +
	geom_point(aes(x=lon_plot[status_plot=="Corrected with alternative method(s)"], y=lat_plot[status_plot=="Corrected with alternative method(s)"]), color=my_palette[3], size=0.7) +
	geom_point(aes(x=lon_plot[status_plot=="Not corrected"], y=lat_plot[status_plot=="Not corrected"]), color=my_palette[2], size=0.7) +
	annotate("rect", xmin=lon_range[1], xmax=lon_range[2], ymin=lat_range[1], ymax = lat_range[2], color="black", fill="black", alpha=0.5) +
	scale_x_continuous(breaks=seq(-180,180,30)) +
	scale_y_continuous(breaks=seq(-90,90,30)) +
	theme_bw() +
	labs(color = "Correction status", x="Longitude", y="Latitude") 
g3.5_all = ggplot() +
	geom_polygon(data=map_world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot, y=lat_plot), color=my_palette[4], size=0.7) +
	annotate("rect", xmin=lon_range[1], xmax=lon_range[2], ymin=lat_range[1], ymax = lat_range[2], color="black", fill="black", alpha=0.5) +
	scale_x_continuous(breaks=seq(-180,180,30)) +
	scale_y_continuous(breaks=seq(-90,90,30)) +
	theme_bw() +
	labs(color = "Correction status", x="Longitude", y="Latitude") 
g4 = ggplot() +
	geom_polygon(data=map_world, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot[status_plot=="Corrected with alternative method(s)"], y=lat_plot[status_plot=="Corrected with alternative method(s)"]), color=my_palette[3], size=0.7) +
	geom_point(aes(x=lon_plot[status_plot=="Not corrected"], y=lat_plot[status_plot=="Not corrected"]), color=my_palette[2], size=0.7) +
	annotate("rect", xmin=lon_range[1], xmax=lon_range[2], ymin=lat_range[1], ymax = lat_range[2], color="black", fill="black", alpha=0.5) +
	scale_x_continuous(breaks=seq(-180,180,30)) +
	scale_y_continuous(breaks=seq(-90,90,30)) +
	theme_bw() +
	labs(color = "Correction status", x="Longitude", y="Latitude") 
g5 = ggplot() +
	geom_polygon(data=map_mb, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Corrected"], y=lat_plot_mb[status_plot_mb=="Corrected"]), color=my_palette[1], size=0.7) +
	theme_bw() +
	theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
g5.5 = ggplot() +
	geom_polygon(data=map_mb, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Corrected"], y=lat_plot_mb[status_plot_mb=="Corrected"]), color=my_palette[1], size=0.7) +
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Corrected with alternative method(s)"], y=lat_plot_mb[status_plot_mb=="Corrected with alternative method(s)"]), color=my_palette[3], size=0.7) +
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Not corrected"], y=lat_plot_mb[status_plot_mb=="Not corrected"]), color=my_palette[2], size=0.7) +
	theme_bw() +
	theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
g5.5_all = ggplot() +
	geom_polygon(data=map_mb, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot_mb, y=lat_plot_mb), color=my_palette[4], size=0.7) +
	theme_bw() +
	theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
g6 = ggplot() +
	geom_polygon(data=map_mb, aes(x=long, y=lat, group=group), fill="#dddddd")+
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Corrected with alternative method(s)"], y=lat_plot_mb[status_plot_mb=="Corrected with alternative method(s)"]), color=my_palette[3], size=0.7) +
	geom_point(aes(x=lon_plot_mb[status_plot_mb=="Not corrected"], y=lat_plot_mb[status_plot_mb=="Not corrected"]), color=my_palette[2], size=0.7) +
	theme_bw() +
	theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())


#png(filename="~/Documents/radiometry_QC/profile_map.png", width=700, height=400)
#plot(g1)
#dev.off()
#png(filename="~/Documents/radiometry_QC/profile_map_medblack.png", width=200, height=150)
#plot(g2)
#dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_corrected_profile_map.png", width=700, height=400)
plot(g3)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_not_corrected_profile_map.png", width=700, height=400)
plot(g4)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_corrected_profile_map_medblack.png", width=220, height=120)
plot(g5)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_not_corrected_profile_map_medblack.png", width=220, height=120)
plot(g6)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_all_profile_map.png", width=700, height=400)
plot(g3.5)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_all_profile_map_medblack.png", width=220, height=120)
plot(g5.5)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_all_profile_map_same_color.png", width=700, height=400)
plot(g3.5_all)
dev.off()
png(filename="~/Documents/radiometry_QC/good_floats_all_profile_map_medblack_same_color.png", width=220, height=120)
plot(g5.5_all)
dev.off()
