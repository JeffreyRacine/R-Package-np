# Shape of a serialized control from 0.80-2. Fresh controls have no calibration.
historical_boot_control <- function(calibration, ...) {
  control <- np_boot_control(...)
  control$calibration <- calibration
  control
}
