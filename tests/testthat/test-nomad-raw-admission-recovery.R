test_that("single native evaluation status distinguishes raw validity from penalties", {
  valid <- list(eval.history=1,invalid.history=0,fval=.5)
  expect_true(.np_nn_single_eval_admissible(valid))
  bad <- valid; bad$invalid.history <- 1
  expect_false(.np_nn_single_eval_admissible(bad))
  for(value in c(NA_real_,Inf,-Inf,.Machine$double.xmax,-.Machine$double.xmax)) {
    bad <- valid; bad$fval <- value
    expect_false(.np_nn_single_eval_admissible(bad))
  }
  for(history in list(numeric(),c(1,1),0,2,NA_real_)) {
    bad <- valid; bad$eval.history <- history
    expect_error(.np_nn_single_eval_admissible(bad),"single native evaluation")
  }
  for(history in list(numeric(),c(0,0),-1,2,NA_real_)) {
    bad <- valid; bad$invalid.history <- history
    expect_error(.np_nn_single_eval_admissible(bad),"single native evaluation")
  }
})

test_that("shared callback admission preserves exploration and permits one recovery", {
  skip_on_cran()
  skip_if_not_installed("crs",minimum_version="0.15.46")
  if(!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  driver <- function(evaluator,recover=NULL,nmulti=1L,remin=FALSE) {
    .np_nomad_search(engine="nomad",baseline_record=NULL,start_degree=1L,
      x0=0,bbin=1L,lb=0,ub=4,eval_fun=evaluator,
      build_payload=function(point,best_record,solution,interrupted)
        list(payload=list(point=point,objective=best_record$objective)),
      native.r.bridge=TRUE,recover_start=recover,nmulti=nmulti,remin=remin,
      nomad.opts=list(MAX_BB_EVAL=20L))
  }
  plain <- function(point) list(objective=(point[1L]-3)^2+1,degree=1L,num.feval=1L)
  explicit <- function(point) c(plain(point),list(admissible=TRUE))
  a <- driver(plain); b <- driver(explicit)
  for(field in c("best_point","n.unique","n.visits","best.restart","restart.starts",
                 "best_payload","native.diagnostics"))
    expect_identical(a[[field]],b[[field]])
  expect_false("admissible" %in% names(b$trace))
  expect_false("admissible" %in% names(b$best))
  phase <- FALSE; recoveries <- 0L; transcript <- numeric()
  evaluator <- function(point) {
    value <- if(phase) (point[1L]-3)^2+1 else 10
    transcript <<- c(transcript,value)
    list(objective=value,degree=1L,num.feval=1L,admissible=phase)
  }
  recover <- function(start) {
    expect_identical(start,0)
    recoveries <<- recoveries+1L
    phase <<- TRUE
    list(found=TRUE,point=3,evaluations=1L,objective=1)
  }
  recovered <- driver(evaluator,recover)
  expect_identical(recoveries,1L)
  expect_length(recovered$restart.results,2L)
  expect_true(recovered$restart.results[[2L]]$recovery)
  expect_identical(recovered$restart.results[[2L]]$start,3)
  expect_identical(recovered$best_point,3)
  expect_equal(recovered$best$objective,1)
  expect_identical(transcript[1L],10)
  expect_identical(recovered$restart.results[[1L]]$objective,10)
  phase <- FALSE; recoveries <- 0L
  repeated <- driver(evaluator,recover,nmulti=2L,remin=TRUE)
  expect_identical(recoveries,1L)
  expect_length(repeated$restart.results,4L)
  expect_true(repeated$restart.results[[3L]]$recovery)
  expect_true(repeated$restart.results[[4L]]$remin)
  expect_identical(repeated$nomad.remin.index,4L)
  expect_identical(repeated$best_point,3)
  never <- function(point) list(objective=10,degree=1L,admissible=FALSE)
  expect_error(driver(never),"failed to obtain any admissible fitted model")
  for(flag in list(1,NA,c(TRUE,FALSE))) {
    wrong <- function(point) list(objective=1,degree=1L,admissible=flag)
    expect_error(driver(wrong),"admissible.*TRUE or FALSE")
  }
  expect_error(driver(never,function(start)list(found=NA)),"scalar logical")
  expect_error(driver(never,function(start)list(found=TRUE,point=4.5)),"invalid start")
  expect_error(driver(never,function(start)list(found=TRUE,point=2.5)),"invalid start")
})
