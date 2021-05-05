library(stringr)
library(ggplot2)
library(latex2exp)

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

print(length(which(-df$mean/100 >= 0.08 & -df$mean/100 <=0.12))/length(df$mean))
print(length(which(-df2$mean/100 >= 0.08 & -df2$mean/100 <=0.12))/length(df2$mean))

p <- ggplot() +
    geom_histogram( aes(x = -df2$mean/100, y=0.005*length(df2$mean)*..density..), fill="#69b3a2" , color="#FFFFFF", boundary=T, binwidth=0.005  ) +
	scale_x_continuous(limits = c(0.05, 0.13), breaks=seq(0.06, 0.14, 0.02)) +
	geom_vline(xintercept=0.08, linetype="dashed") +
	geom_vline(xintercept=0.12, linetype="dashed") +
    #labs(x="- dP/dt (decibar/s)", y="Number of profiles") +
    #labs(x=TeX("$-\\frac{dP}{dt}$ ($dbar\\cdot s^{-1}$)"), y="Number of profiles") +
    labs(x=TeX("Ascending speed (dbar s^{-1})"), y="Number of profiles") +
	theme_bw()

png(filename="~/Documents/radiometry_QC/pretty_plots/asc_speeds.png", width=300, height=300)
p
dev.off()
