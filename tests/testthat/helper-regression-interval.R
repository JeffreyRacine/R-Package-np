# The fixtures parameterize internal providers by their calibration names.
regression_test_interval <- function(calibration) {
  if (is.null(calibration)) return(NULL)
  switch(calibration, absolute="symmetric-sd", "separate-tails"="asymmetric-sd",
    "tail-envelope"="asymmetric-sd-envelope", calibration)
}
