test_that("npudens all-categorical tree profile bandwidth CV matches dense CV", {
  if (exists("spawn_mpi_slaves", mode = "function"))
    spawn_mpi_slaves(1L)

  old.tree <- getOption("np.tree")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.messages = old.messages)
    if (exists("close_mpi_slaves", mode = "function"))
      close_mpi_slaves()
  }, add = TRUE)
  options(np.messages = FALSE)

  run_case <- function(kind, bwmethod, ukertype, okertype, seed) {
    set.seed(seed)
    n <- 300L

    if (identical(kind, "unordered")) {
      dat <- data.frame(
        u1 = factor(sample(letters[1:4], n, TRUE)),
        u2 = factor(sample(LETTERS[1:3], n, TRUE))
      )
      form <- ~ u1 + u2
    } else if (identical(kind, "ordered")) {
      dat <- data.frame(
        o1 = ordered(sample(1:6, n, TRUE)),
        o2 = ordered(sample(1:5, n, TRUE))
      )
      form <- ~ o1 + o2
    } else {
      dat <- data.frame(
        u1 = factor(sample(letters[1:4], n, TRUE)),
        o1 = ordered(sample(1:6, n, TRUE))
      )
      form <- ~ u1 + o1
    }

    options(np.tree = FALSE)
    bw.dense <- npudensbw(
      form,
      data = dat,
      nmulti = 1,
      bwmethod = bwmethod,
      ukertype = ukertype,
      okertype = okertype
    )

    options(np.tree = TRUE)
    bw.profile <- npudensbw(
      form,
      data = dat,
      nmulti = 1,
      bwmethod = bwmethod,
      ukertype = ukertype,
      okertype = okertype
    )

    expect_equal(
      bw.profile$bw,
      bw.dense$bw,
      tolerance = 1e-8,
      info = paste(kind, bwmethod, ukertype, okertype, "bw")
    )
    expect_equal(
      bw.profile$fval,
      bw.dense$fval,
      tolerance = 1e-8,
      info = paste(kind, bwmethod, ukertype, okertype, "fval")
    )
    expect_equal(
      bw.profile$num.feval,
      bw.dense$num.feval,
      info = paste(kind, bwmethod, ukertype, okertype, "num.feval")
    )
  }

  run_case("unordered", "cv.ml", "aitchisonaitken", "wangvanryzin", 2001L)
  run_case("unordered", "cv.ls", "liracine", "wangvanryzin", 2002L)
  run_case("ordered", "cv.ml", "aitchisonaitken", "racineliyan", 2003L)
  run_case("ordered", "cv.ls", "aitchisonaitken", "liracine", 2004L)
  run_case("mixed", "cv.ml", "liracine", "racineliyan", 2005L)
  run_case("mixed", "cv.ls", "liracine", "liracine", 2006L)
})
