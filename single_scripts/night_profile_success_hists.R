library(ggplot2)

n_nights_success = c(4, 3, 8, 8, 1, 1, 32, 1, 10, 1, 1, 5, 2, 4, 1, 1, 3, 12, 2, 7, 1, 2, 16, 1, 2, 3, 3, 17, 5, 6, 1, 3, 5, 1, 1, 1, 1, 5, 14, 4, 14, 1, 2, 2, 5, 14)

n_nights_fail = c(3, 1, 2, 3, 1, 1, 2, 3, 1)

df = data.frame(state = c(rep("Method succeeded", length(n_nights_success)), rep("Method failed", length(n_nights_fail))),
				value = c(n_nights_success, n_nights_fail))

p = ggplot(df, aes(x=value, fill=state)) +
	#geom_bar(color="#e9ecef", alpha=0.6, position='identity') +
	geom_bar(color="#e9ecef") +
	#geom_bar(color="#e9ecef", alpha=0.6, position=position_dodge2(preserve="single")) +
	theme_bw() +
	scale_y_continuous(limits=c(0,20), breaks=seq(0,20,5), minor_breaks=seq(0,20)) +
	scale_x_continuous(limits=c(0,35), breaks=seq(0,35,5)) +
	xlab("Number of available night profiles per float") +
	ylab("Number of floats") +
	labs(fill="") +
	theme(legend.position = c(.95, .95),
		  legend.justification = c("right", "top"),
		  legend.box.just = "right",
		  legend.margin = margin(6,6,6,6))#,
		  #legend.box.background = element_rect(fill="white", colour="black"))

png(filename="night_profiles_hist.png", width=400, height=200)
p
dev.off()
