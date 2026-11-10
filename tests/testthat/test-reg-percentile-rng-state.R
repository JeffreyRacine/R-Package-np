test_that("regression percentile tasks replay the coordinator generator", {
  skip_if(getRversion() < "3.6.0")
  kind <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({do.call(RNGkind, as.list(kind)); if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  state <- .np_reg_percentile_rng_state(91)
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = "."),
    header = state[1L], length = length(state), kind = RNGkind())
  expected <- .np_reg_percentile_with_seed(91, runif(37))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(17)
  worker <- .Random.seed
  expect_identical(.np_reg_percentile_with_state(state, contract, runif(37)), expected)
  expect_identical(.Random.seed, worker)
  expect_error(.np_reg_percentile_with_state(NULL, contract, stop("entered")), "invalid regression percentile")
  expect_error(.np_reg_percentile_with_state(state, contract, stop("primary")), "primary")
  expect_identical(.Random.seed, worker)
  state[20L] <- NA_integer_
  expect_silent(.np_reg_percentile_rng_validate(state, contract))
})


test_that("task RNG cleanup survives throwing warning handlers", {
  skip_if(getRversion() < "3.6.0")
  kind <- RNGkind(); had <- exists('.Random.seed', .GlobalEnv, inherits = FALSE)
  old <- if (had) get('.Random.seed', .GlobalEnv) else NULL
  on.exit({suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (had) assign('.Random.seed', old, .GlobalEnv) else
      if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv)})
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  state <- .np_reg_percentile_rng_state(91)
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = '.'),
    header = state[1L], length = length(state), kind = RNGkind())
  bad <- state; bad[1L] <- 99999L; invalid <- contract; invalid$header <- bad[1L]
  expect_error(.np_reg_percentile_rng_validate(bad, invalid), 'invalid regression percentile')
  invalid <- contract; invalid$kind[1L] <- 'Super-Duper'
  expect_error(.np_reg_percentile_rng_validate(state, invalid), 'invalid regression percentile')
  for (absent in c(FALSE, TRUE)) {
    suppressWarnings(RNGkind("L'Ecuyer-CMRG", 'Inversion', 'Rounding'))
    set.seed(42)
    expected.header <- .Random.seed[1L]
    if (absent) rm('.Random.seed', envir = .GlobalEnv)
    entry <- if (!absent) .Random.seed else NULL
    expected.kind <- kind
    expected.kind[1:3] <- c("L'Ecuyer-CMRG", 'Inversion', 'Rounding')
    err <- tryCatch(withCallingHandlers(
      .np_reg_percentile_with_state(state, contract, {runif(3); warning('body')}),
      warning = function(w) stop('handler primary')), error = identity)
    expect_match(conditionMessage(err), 'handler primary')
    expect_identical(if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) .Random.seed else NULL, entry)
    if (!absent) expect_identical(RNGkind(), expected.kind)
    else expect_identical(.np_rng_scope_capture()$seed[1L], expected.header)
  }
})


test_that("task RNG validation respects the runtime kind layout", {
  # R before 3.6 has two kind fields and the Rounding sampler is implicit.
  validate <- .np_reg_percentile_rng_validate
  environment(validate) <- list2env(list(.np_rng_kind_count = function() 2L),
    parent = environment(validate))
  state <- c(403L, rep.int(1L, 625L))
  contract <- list(runtime = paste(R.version$major, R.version$minor, sep = '.'),
    header = state[1L], length = length(state), kind = c('Mersenne-Twister', 'Inversion'))
  expect_silent(validate(state, contract))
  contract$kind <- c(contract$kind, 'Rounding')
  expect_error(validate(state, contract), 'invalid regression percentile')
})


test_that("task contracts distinguish old and fourth-kind layouts", {
  validate <- .np_reg_percentile_rng_validate
  environment(validate) <- list2env(list(.np_rng_kind_count=function() 4L),parent=environment(validate))
  state <- c(110403L,rep.int(1L,625L))
  contract <- list(runtime=paste(R.version$major,R.version$minor,sep='.'),
    header=state[1L],length=length(state),kind=c('Mersenne-Twister','Inversion','Rejection','BTPE'))
  expect_silent(validate(state,contract))
  contract$kind[4L]<-'unknown'
  expect_error(validate(state,contract),'invalid regression percentile')
  contract$kind[4L]<-'Buggy BTPE'
  expect_error(validate(state,contract),'invalid regression percentile')
  state[1L]<-10403L;contract$header<-state[1L]
  expect_silent(validate(state,contract))
  contract$kind<-contract$kind[1:3]
  expect_error(validate(state,contract),'invalid regression percentile')
})

test_that("runtime built-in kinds preserve task draws and caller continuation", {
  skip_if(getRversion() < '3.6.0')
  old.kind<-RNGkind(); had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  old<-if(had) .Random.seed else NULL
  on.exit({suppressWarnings(do.call(RNGkind,as.list(old.kind)))
    if(had) assign('.Random.seed',old,.GlobalEnv) else
      if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)},add=TRUE)
  uniforms<-c('Wichmann-Hill','Marsaglia-Multicarry','Super-Duper','Mersenne-Twister',
    'Knuth-TAOCP','Knuth-TAOCP-2002',"L'Ecuyer-CMRG")
  binoms<-if(.np_rng_kind_count()==4L) c('Buggy BTPE','BTPE') else NA_character_
  draw<-function() list(runif(9),rnorm(7),sample.int(100,11),rbinom(13,100,.3))
  for(u in uniforms) for(n in c('Buggy Kinderman-Ramage','Ahrens-Dieter','Inversion','Kinderman-Ramage'))
    for(s in c('Rounding','Rejection')) for(b in binoms) {
      kinds<-list(kind=u,normal.kind=n,sample.kind=s)
      if(!is.na(b)) kinds$binom.kind<-b
      suppressWarnings(do.call(RNGkind,kinds));set.seed(73)
      entry<-.Random.seed; expected<-draw(); continuation<-.Random.seed
      assign('.Random.seed',entry,.GlobalEnv)
      contract<-list(runtime=paste(R.version$major,R.version$minor,sep='.'),
        header=entry[1L],length=length(entry),kind=RNGkind())
      expect_identical(.np_reg_percentile_with_state(entry,contract,draw()),expected)
      expect_identical(.Random.seed,entry)
      expect_error(.np_reg_percentile_with_state(entry,contract,{draw();stop('body error')}),'body error')
      expect_identical(.Random.seed,entry)
      expect_identical(draw(),expected)
      expect_identical(.Random.seed,continuation)
    }
})
