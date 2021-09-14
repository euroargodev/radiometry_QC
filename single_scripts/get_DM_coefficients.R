require(ncdf4)
require(stringr)
require(ggplot2)
require(gridExtra)

source("~/Documents/radiometry/helpers.R")

WMO_list = read.table("~/Documents/radiometry_QC/floats_corrected_main_radiometry.txt")$V1

df = NULL

df$WMO = rep(WMO_list, each=4)

PARAM_NAMES = c("DOWN_IRRADIANCE380", "DOWN_IRRADIANCE412", 
                "DOWN_IRRADIANCE490", "DOWNWELLING_PAR")

irr_labels = c("E[d](380)", "E[d](412)", "E[d](490)", "PAR")
names(irr_labels) = PARAM_NAMES

df$PARAM = rep(PARAM_NAMES, length(WMO_list))

n = length(df$WMO)
df$A = rep(NA, n)
df$B = rep(NA, n)
df$C = rep(NA, n)
df$Q = rep(NA, n)

for (i in 1:n) {
  filename = paste0("/mnt/c/DATA/ftp.ifremer.fr/ifremer/argo/dac/coriolis/", 
                    df$WMO[i], "/profiles/BD", df$WMO[i], "_001.nc")
  fnc = nc_open(filename)
  coef = ncvar_get(fnc, "SCIENTIFIC_CALIB_COEFFICIENT")  
  param = ncvar_get(fnc, "PARAMETER")
  nc_close(fnc)
  id = which(param == str_pad(df$PARAM[i], 64, side="right"))
  calib = unlist(str_split(coef[id], ","))
  for (letter in c("A", "B", "C", "Q")) {
    if (any(grepl(letter, calib))) {
      value = unlist(str_split(calib[grepl(letter, calib)], "="))[2]
      df[[letter]][i] = as.numeric(value)
    } else {
      df[[letter]][i] = 0
    }
  }
}

df = data.frame(df)

df2 = df[which(df$Q==0),]

gA = ggplot(data=df2, mapping=aes(x=A, group=PARAM)) +
  geom_histogram() +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE),
                     guide = guide_axis(check.overlap=TRUE)) +
  theme_bw() +
  #theme(axis.text=element_text(size=9)) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~PARAM, ncol=1, scales="free", 
             labeller=labeller(PARAM=as_labeller(irr_labels, default=label_parsed)))

gB = ggplot(data=df2, mapping=aes(x=B, group=PARAM)) +
  geom_histogram() +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE),
                     guide = guide_axis(check.overlap=TRUE)) +
  theme_bw() +
  #theme(axis.text=element_text(size=9)) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~PARAM, ncol=1, scales="free", 
             labeller=labeller(PARAM=as_labeller(irr_labels, default=label_parsed)))

gC = ggplot(data=df2, mapping=aes(x=C, group=PARAM)) +
  geom_histogram() +
  scale_y_continuous(breaks=c(0,2,4,6,8)) +
  scale_x_continuous(labels = function(x) format(x, scientific=TRUE),
                     guide = guide_axis(check.overlap=TRUE)) +
  theme_bw() +
  #theme(axis.text=element_text(size=9)) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~PARAM, ncol=1, scales="free", 
             labeller=labeller(PARAM=as_labeller(irr_labels, default=label_parsed)))

png(filename="~/Documents/radiometry_QC/DM_coefficients_hists.png", width=800, height=800)
grid.arrange(gA, gB, gC, nrow=1, left="count")
dev.off()
#gg1 = arrangeGrob(gA, gB, gC, nrow=1, left="count")
#ggsave(filename="~/Documents/radiometry_QC/DM_coefficients_hists.tiff", plot=gg1, device="tiff", width=180, height=180, units="mm", dpi=300)
