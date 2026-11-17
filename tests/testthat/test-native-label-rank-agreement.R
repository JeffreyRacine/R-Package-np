test_that("native names carry coordinator text and omit identity-owned roles", {
  ns <- asNamespace("npRmpi")
  plan <- get(".npRmpi_native_name_plan", ns)
  apply <- get(".npRmpi_native_names_apply", ns)
  unpack <- get(".npRmpi_native_label_unpack", environment(plan))
  old <- options(scipen = 999)
  on.exit(options(old), add = TRUE)
  make <- function(xdat, ydat, zdat, identity.roles = character()) {
    plan(sys.function(), environment(), parent.frame(),
         list(xdat = xdat, ydat = ydat, zdat = zdat), identity.roles)
  }
  x <- seq_len(10L)
  p <- make(x * 1e5, x * 1e-4, x * 1e6)
  expect_identical(p$xnames, paste(deparse(quote(x * 1e5)), collapse = ""))
  expect_identical(p$ynames, deparse(quote(x * 1e-4)))
  expect_identical(p$znames, paste(deparse(quote(x * 1e6)), collapse = ""))
  expect_true(all(vapply(p, is.character, logical(1L))))
  b <- structure(list(xnames = ".__npRmpi_autod_x", ynames = ".__npRmpi_autod_y",
                      znames = ".__npRmpi_autod_z"), class = "rbandwidth")
  coordinator <- apply(b, p)
  options(scipen = 0)
  worker <- apply(b, unserialize(serialize(p, NULL)))
  expect_identical(worker, coordinator)
  expect_identical(worker$xnames, "x * 100000")
  expect_identical(worker$ynames, "x * 0.0001")
  expect_identical(make(x, x, x, c("xdat", "ydat", "zdat")), list())
  expect_identical(names(make(x, x, x, "tydat")), c("ynames", "xnames", "znames"))
  expect_identical(names(make(x, x, x, "ydat")), c("xnames", "znames"))
  b$xnames <- "already supplied"
  expect_identical(apply(b, p)$xnames, "already supplied")
  long <- as.call(c(list(as.name("c")), as.list(seq_len(100L))))
  expected <- deparse(long)
  expect_gt(length(expected), 1L)
  owner <- list2env(list(xdat = quote(x), ydat = long), parent = baseenv())
  actual <- plan(function(xdat, ydat) NULL, owner, baseenv(), list(xdat = x, ydat = x))
  expect_identical(unpack(actual$ynames), expected)
  for (name in c("dat", "tdat")) {
    definition <- function(dat) NULL
    names(formals(definition)) <- name
    owner <- list2env(setNames(list(quote(x * 1e5)), name), parent = baseenv())
    result <- plan(definition, owner, baseenv(), setNames(list(x), name))
    expect_identical(result$xnames, paste(deparse(quote(x * 1e5)), collapse = ""))
  }
  for (label in c("caf\u00e9", enc2native("caf\u00e9"), iconv("caf\u00e9", to="latin1"))) {
    p <- list(ynames = label)
    expect_identical(apply(b,p)$ynames, apply(b,unserialize(serialize(p,NULL)))$ynames)
  }
})

test_that("large exact native labels use only smaller validated encodings", {
  plan <- getFromNamespace(".npRmpi_native_name_plan", "npRmpi")
  pack <- get(".npRmpi_native_label_pack", environment(plan))
  unpack <- get(".npRmpi_native_label_unpack", environment(plan))
  offset <- length(serialize("", NULL))
  for (size in c(511L, 512L, 513L)) {
    text <- strrep("a", size - offset)
    expect_identical(length(serialize(text, NULL)), size)
    encoded <- pack(text)
    expect_identical(unpack(encoded), text)
    expect_identical(is.character(encoded), size < 512L)
    expect_lte(length(serialize(encoded, NULL)), size)
  }
  label <- c(strrep("x", 2000L), strrep("y", 2000L))
  encoded <- pack(label)
  expect_type(encoded, "list")
  expect_identical(unpack(unserialize(serialize(encoded, NULL))), label)
  expect_error(unpack(list(codec = "unknown", payload = raw())), "invalid native label metadata")
  expect_error(unpack(list(codec = "bzip2", payload = 1)), "invalid native label metadata")
  bad <- list(codec = "bzip2", payload = memCompress(serialize(1, NULL), "bzip2"))
  expect_error(unpack(bad), "invalid native label text")
  bad$payload <- memCompress(serialize(setNames("a", "name"), NULL), "bzip2")
  expect_error(unpack(bad), "invalid native label text")
  expect_error(unpack(NA_character_), "invalid native label text")
  expect_error(unpack(character()), "invalid native label text")
  bad$payload <- memCompress(serialize(list(a = 1), NULL), "bzip2")
  expect_error(unpack(bad), "invalid native label text")
  bad$payload <- memCompress(as.raw(1:4), "bzip2")
  expect_error(unpack(bad))
  bad$payload <- encoded$payload[seq_len(length(encoded$payload) - 8L)]
  expect_error(suppressWarnings(unpack(bad)))
})
