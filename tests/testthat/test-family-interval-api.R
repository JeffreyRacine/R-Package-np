test_that("explicit family intervals preserve distinct fitted providers", {
  base <- list(errors = "bootstrap", center = "bias-corrected", band = "all",
    bootstrap = "inid", interval = "asymmetric-sd")
  translated <- .np_plot_family_interval_request(base, "npudens")
  expect_identical(translated$plot.errors.boot.calibration, "asymmetric-sd")
  expect_null(translated$interval)
  expect_identical(translated$.np.interval.request$interval, "asymmetric-sd")
  both <- base
  both$plot.errors.boot.calibration <- "absolute"
  expect_error(.np_plot_family_interval_request(both, "npudens"), "replace")
  scoped <- base
  scoped$band <- NULL
  scoped$bootstrap <- "smoothed-q6"
  expect_error(.np_plot_family_interval_request(scoped, "npudens"), "explicit band")
  expect_error(.np_plot_check_public_names("plot.errors.boot.calibration"), "use interval =", fixed = TRUE)
  marginal <- list(errors = "bootstrap", interval = "normal-sd")
  expect_identical(.np_plot_family_interval_request(marginal, "npudens")$band, "pmzsd")
  ordinary <- list(errors = "bootstrap", interval = "percentile", band = "pointwise")
  expect_identical(.np_plot_family_interval_request(ordinary, "npudens")$band, "pointwise")
  old <- np_boot_control(blocklen = 3L)
  old$calibration <- "default"
  expect_error(.np_plot_family_interval_request(list(boot.control = old), "npudist"), "rebuild")
  expect_error(np_boot_control("exact", "mammen", 3L, "absolute"), "move calibration")
})

test_that("component metadata reflects the provider and returned rows", {
  malformed <- list(schema_version = 2L, family = c("npudens", "npudist"),
                    request = list(interval = NULL))
  expect_error(.np_plot_data_validate(malformed), "invalid np_plot_data: unknown schema or family",
               fixed = TRUE)
  rows <- data.frame(band = c("pointwise", "simultaneous"), status = c("finite", "undefined"))
  describe <- function(actual, family = "npudens") .np_plot_family_components(
    list(construction = "density-joint-scale-envelope-S-v1", calibration = actual),
    rows, list(errors = "bootstrap"), family)$components
  expect_identical(describe("separate-tails")$pointwise$interval, "separate-tails")
  expect_identical(describe("separate-tails", "npudist")$pointwise$interval, "asymmetric-sd")
  expect_identical(describe("percentile-tails")$pointwise$interval, "percentile-tails")
  expect_identical(describe("asymmetric-sd")$pointwise$interval, "asymmetric-sd")
  expect_identical(describe("asymmetric-sd")$simultaneous$interval, NA_character_)
})

test_that("smoothed beta-binomial accepts the ordinary bootstrap control default", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  local <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(local), add = TRUE)
  x <- data.frame(x = seq(.04, .96, length.out = 36L)^1.3)
  b <- npudensbw(dat = x, bws = .14, bandwidth.compute = FALSE,
    ckertype = "epanechnikov", ckerorder = 4L, ckerbound = "fixed", ckerlb = 0, ckerub = 1)
  p <- list(x = b, xdat = x, center = "bias-corrected", errors = "bootstrap",
    bootstrap = "smoothed-beta-binomial", B = 19L, random.seed = 49L,
    interval = "symmetric-sd", band = "pointwise", neval = 4L, output = "data")
  implicit <- suppressWarnings(do.call(plot, p))
  p$boot.control <- np_boot_control()
  explicit <- suppressWarnings(do.call(plot, p))
  expect_identical(implicit$panels, explicit$panels)
  expect_identical(implicit$inference[[1L]]$resampling_law, "smoothed-beta-binomial")
  expect_identical(implicit$inference[[1L]]$components$pointwise$interval, "symmetric-sd")
})
