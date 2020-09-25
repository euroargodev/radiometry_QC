library(stringr)
library(ggplot2)

# CYCLE_NUMBER, median, mean, max, min, sd

file_list = system2("ls", args=c("~/Documents/radiometry_QC/asc_speeds/*"), stdout = T)

WMO_list = rep(NA, length(file_list))
for (i in 1:length(file_list)) {
    separated = unlist(str_split(file_list[i], "/"))
    WMO_list[i] = str_sub(separated[length(separated)], 1, 7)
}

df = data.frame("WMO"=character(0), "CYCLE_NUMBER"=numeric(0), "median"=numeric(0), "mean"=numeric(0), "max"=numeric(0), "min"=numeric(0), "sd"=numeric(0))

for (i in 1:length(file_list)) {
    open_tab = read.table(file_list[i], col.names = c("CYCLE_NUMBER", "median", "mean", "max", "min", "sd"))
    open_tab$WMO = rep(WMO_list[i], dim(open_tab)[1])
    
    df = rbind(df, open_tab)
}

df$mean[!is.finite(df$mean)] = NA
df$max[!is.finite(df$max)] = NA

#df[which(df$mean>10000),]$max = NA
#df[which(df$mean>10000),]$min = NA
#df[which(df$mean>10000),]$sd = NA
#df[which(df$mean>10000),]$mean = NA

for (param in c("median", "mean", "max", "min", "sd")) {
    param_axis = df[[param]]
    quant = quantile(param_axis, na.rm=T, probs = c(1/4, 1/2, 3/4))
    
    outlier = which(param_axis < quant[2] - 7*(quant[3] - quant[1]) | param_axis > quant[2] + 7*(quant[3] - quant[1]))
    df[[param]][outlier] = NA
}

df2 = na.omit(df)

p <- ggplot() +
    geom_histogram( aes(x = df2$mean, y=0.5*length(df2$mean)*..density..), fill="#69b3a2" , color="#FFFFFF", boundary=T, binwidth=0.5  ) +
    #geom_label( aes(x=250, y=6, label="corrected"), color="#69b3a2") +
    labs(x="dP/dt (decibar/s)", y="Number of profiles")
