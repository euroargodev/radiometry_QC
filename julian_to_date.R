require(date)

solar_angle_test <- function (julian_day, latitude, longitude) {

  # get the number of days without the hour
  nb_day = floor(julian_day)

  # Then calculate the number of hours 
  nb_hour = julian_day-nb_day

  # Round it 
  hour = round(24*nb_hour)

  # get the month and the day from the julian day
  jd = date.mdy(julian_day)

  # the day
  day = jd$day

  # the month
  month = jd$month

  # calculate the solar angle
  asol = possol.vec(month, day, hour, longitude, latitude)

  asol_zen = asol$zen

  solar_test_bool = (asol_zen >= 2.)

  return(solar_test_bool)

}
