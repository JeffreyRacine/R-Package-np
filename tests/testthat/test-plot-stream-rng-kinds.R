test_that("isolated plot request tails recognize only supported runtime kinds", {
  kinds <- c("Mersenne-Twister", "Inversion", "Rejection")
  absent <- !exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if (!absent) seed <- .Random.seed
  binoms <- if (.np_rng_kind_count() == 4L) c("Buggy BTPE", "BTPE") else character()
  for (tail in c(list(character()), lapply(binoms, identity))) {
    k <- c(kinds,tail)
    expect_identical(.np_udens_beta_random_request("multinomial",10,7,42,k)$kind,k)
    k[2] <- "Box-Muller"
    expect_identical(.np_plot_normal_request(42,k,7)$kind,k)
  }
  for (k in list(c(kinds,"unknown"),c(kinds,"BTPE","future"),
                 c(kinds[1:2],"unknown"),c(kinds,NA_character_))) {
    expect_error(.np_udens_beta_random_request("normal",10,7,42,k),"invalid beta")
    k[2] <- "Box-Muller"
    expect_error(.np_plot_normal_request(42,k,7),"supported RNG kinds")
  }
  if (.np_rng_kind_count() < 4L) {
    expect_error(.np_udens_beta_random_request("normal",10,7,42,c(kinds,"BTPE")),"invalid beta")
    expect_error(.np_plot_normal_request(42,c(kinds[1],"Box-Muller",kinds[3],"BTPE"),7),"supported RNG kinds")
  }
  expect_identical(!exists(".Random.seed",.GlobalEnv,inherits=FALSE),absent)
  if (!absent) expect_identical(.Random.seed,seed)
})

test_that("isolated plot samplers preserve complete runtime kinds and caller continuation", {
  old <- RNGkind();present <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  saved <- if (present) .Random.seed else NULL
  on.exit({suppressWarnings(do.call(RNGkind,as.list(old)));if(present)assign(".Random.seed",saved,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  binoms <- if (.np_rng_kind_count()==4L)c("Buggy BTPE","BTPE") else NA_character_
  for (binom in binoms) for(normal in c("Inversion","Box-Muller")) {
    k <- list(kind="Mersenne-Twister",normal.kind=normal,sample.kind="Rejection")
    if(!is.na(binom)) k$binom.kind <- binom
    suppressWarnings(do.call(RNGkind,k));kind <- RNGkind()
    for(mode in c("normal","multinomial")) {
      set.seed(42);expected <- if(mode=="normal")matrix(rnorm(70),10) else matrix(as.double(rmultinom(7,10,rep(.1,10))),10)
      set.seed(871);invisible(rnorm(1));entry <- .Random.seed;following <- rnorm(7)
      set.seed(871);invisible(rnorm(1))
      store <- .np_udens_beta_store()
      stream <- tryCatch(.np_udens_beta_random_stream(store,mode,10,7,42,kind),error=function(e){.np_udens_beta_store_close(store);stop(e)})
      expect_identical(stream$request$kind,kind);expect_identical(stream$receipt$kind,kind)
      expect_identical(.np_udens_beta_random_read(stream,1,7),expected)
      expect_false(stream$receipt$pid==Sys.getpid())
      expect_false(any(c("np","npRmpi","Rmpi")%in%stream$receipt$namespaces))
      expect_identical(.Random.seed,entry);expect_identical(rnorm(7),following)
      .np_udens_beta_store_close(store);expect_false(dir.exists(store$path))
    }
    if(normal=="Box-Muller") {
      set.seed(42);expected <- rnorm(9)
      set.seed(873);invisible(rnorm(1));entry <- .Random.seed;following <- rnorm(7)
      set.seed(873);invisible(rnorm(1))
      stream <- .np_plot_normal_stream(42,kind,9)
      expect_identical(stream$info$kind,kind)
      expect_identical(c(.np_plot_normal_read(stream,1),.np_plot_normal_read(stream,8)),expected)
      expect_false(stream$info$pid==Sys.getpid())
      expect_false(any(c("np","npRmpi","Rmpi")%in%stream$info$namespaces))
      expect_identical(.Random.seed,entry);expect_identical(rnorm(7),following)
      path <- stream$path;.np_plot_normal_close(stream);expect_false(dir.exists(path))
    }
  }
})
