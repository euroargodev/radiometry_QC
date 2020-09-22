library(ggplot2)

drift_range = c(5e-5, 2e-5, 3e-5, 1e-4, 1e-5, 1e-5, 1.5e-4, 4e-5, 2e-5, 7e-5, 2e-5, 1e-4, 2e-4, 2e-5, 8e-5, 5e-5, 1e-5, 2e-5, 5e-5)
T_range = c(1e-4, 4e-5, 5e-5, 5e-6, 3e-4, 4e-5, 3.5e-4, 2e-4, 2e-4, 7e-5, 3e-5, 2e-4, 3e-4, 3e-5, 2e-4, 3e-5, 2e-4, 2e-4, 2e-4)
id = range(drift_range, T_range)
sensor = c(2.5e-5, 2.5e-5)

p = ggplot() +
    geom_line(aes(x=id, y=id), color="red") +
    geom_line(aes(x=id, y=sensor)) +
    geom_line(aes(x=sensor, y=id)) +
    geom_point(aes(x = T_range, y= drift_range)) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(x = "Range of Irr490 correction for temperature", y = "Range of Irr490 correction for drift")
