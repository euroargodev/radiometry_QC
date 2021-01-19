library(ggplot2)
library(ggExtra)
library(parallel)
library(RColorBrewer)


#corr_list = read.csv("~/Documents/radiometry_QC/float_list_temp.t", sep=" ", header=T)
source("~/Documents/radiometry/single_scripts/count_profiles_corrected.R")
#index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)


corr_list = corr_list[which(corr_list$material!="LIVE"),]

pattern = paste0(corr_list$WMO, "_001.nc")
pattern[126] = paste0(corr_list$WMO[126], "_002.nc") #This one does not have profile 001 ?

n_cores = detectCores()

matches = unlist(mcmapply(grep, pattern, mc.cores=n_cores, USE.NAMES=FALSE, MoreArgs=list(x=index_ifremer$file)))

date_list = index_ifremer$date[matches]
date_list = as.Date(as.character(date_list), format="%Y%m%d%H%M%S", tz="UTC")
#date_list = julian(date_list, origin=as.Date("1950-01-01", tz="UTC"))

n_profiles = corr_list$profiles

status = rep(NA, length(date_list))
status[which(corr_list$corrected==1)] = "Corrected"
status[which(corr_list$corrected==0)] = "Not corrected"
status[which(is.na(corr_list$corrected))] = "TODO"


my_palette = brewer.pal(n=3, name="Dark2")[1:2]

g1 = ggplot() +
	geom_point(aes(x=date_list, y=n_profiles, color=status)) +
	scale_color_manual(values=my_palette) +
	scale_x_date(date_labels = "%b %Y") +
	theme_bw() +
	theme(legend.position="none") +
	labs(color = "Correction status", x="Deployment date", y="Number of measured profiles") 

g2 = ggMarginal(g1, groupColour=TRUE, groupFill=TRUE)

png(filename="~/Documents/radiometry_QC/pretty_plots/date_nprof_density.png", width=400, height=400, bg="white")
g2
dev.off()
