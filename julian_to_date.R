solar_angle_test <- function (julian_day, latitude, longitude) {

require(date)
source("./possol.R")

# get the number of days without the hour
nb_jour=floor(julian_day)

# Then calculate the number of hours 
nb_hour=julian_day-nb_jour

# Round it 
hour=round(24*nb_hour)

# get the month and the day from the julian day
jd=date.mdy(julian_day)

# the day
day=jd$day

# the month
month=jd$month

# calculate the solar angle
asol=possol.vec(month,day,hour,longitude,latitude)

asol_zen=asol$zen

if (asol_zen >=2.) 
{
solar_test_bool=TRUE
} else {
solar_test_bool=FALSE
}

return(solar_test_bool)
}
