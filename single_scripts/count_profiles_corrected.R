library(ggplot2)

corr_list = read.csv("~/Documents/radiometry_QC/float_list_temp.t", sep=" ", header=T)
index_ifremer = read.table("~/Documents/radiometry/argo_bio-profile_index.txt", sep=",", header = T)

files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
ident = strsplit(files,"/") #separate the different roots of the files paths
ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
dac = ident[,1] #retrieve the DAC of all profiles as a vector
wod = ident[,2] #retrieve the WMO of all profiles as a vector
prof_id = ident[,4] #retrieve all profiles  name as a vector
variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file



count_profile = rep(NA, length(corr_list$WMO))
    
for (i in 1:length(corr_list$WMO)) {
    n_prof = length(which(grepl("DOWNWELLING_PAR", variables) & wod==corr_list$WMO[i]))
    count_profile[i] = n_prof
}

corr_list$profiles = count_profile

seen = !is.na(corr_list$corrected)

# alive
n_prof_alive = sum(corr_list$profiles[ which(seen & corr_list$material=="LIVE") ])

# aluminium
n_prof_alum = sum(corr_list$profiles[ which(seen & corr_list$material=="ALU") ])

# not corrected 
n_prof_notcorrected = sum(corr_list$profiles[ which(seen & corr_list$material=="PEEK" & corr_list$corrected==0) ])
n_prof_notcorrected_alum = sum(corr_list$profiles[ which(seen & corr_list$material=="ALU" & corr_list$corrected==0) ])

# corrected
n_prof_corrected = sum(corr_list$profiles[ which(seen & corr_list$material=="PEEK" & corr_list$corrected==1) ])
n_prof_corrected_alum = sum(corr_list$profiles[ which(seen & corr_list$material=="ALU" & corr_list$corrected==1) ])

n_prof_alive
n_prof_alum
n_prof_notcorrected
n_prof_corrected

not_corr = corr_list$profiles[ which(seen & corr_list$corrected==0) ]
corr = corr_list$profiles[ which(seen & corr_list$corrected==1) ]

p <- ggplot() +
    geom_histogram( aes(x = corr, y=50*length(corr)*..density..), fill="#69b3a2", color="#FFFFFF", binwidth=50, boundary=T  ) +
    geom_label( aes(x=250, y=6, label="corrected"), color="#69b3a2") +
    geom_histogram( aes(x = not_corr, y=-50*length(not_corr)*..density..), fill= "#404080", color="#FFFFFF", binwidth=50, boundary=T) +
    geom_label( aes(x=250, y=-4, label="not corrected"), color="#404080") +
    labs(x="Number of radiometry profiles in float", y="Number of floats")
