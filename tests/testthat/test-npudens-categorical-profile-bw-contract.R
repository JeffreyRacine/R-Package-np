test_that("npudens all-categorical profile bandwidth CV stays numerically aligned with dense CV", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
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

    options(np.tree = FALSE, np.categorical.compress = FALSE)
    bw.dense <- npudensbw(
      form,
      data = dat,
      nmulti = 1,
      bwmethod = bwmethod,
      ukertype = ukertype,
      okertype = okertype
    )

    options(np.tree = FALSE, np.categorical.compress = TRUE)
    bw.profile <- npudensbw(
      form,
      data = dat,
      nmulti = 1,
      bwmethod = bwmethod,
      ukertype = ukertype,
      okertype = okertype
    )

    expect_true(all(is.finite(bw.profile$bw)),
                info = paste(kind, bwmethod, ukertype, okertype, "bw finite"))
    expect_true(is.finite(bw.profile$fval),
                info = paste(kind, bwmethod, ukertype, okertype, "fval finite"))

    if (identical(bwmethod, "cv.ml")) {
      expect_true(abs(bw.profile$fval - bw.dense$fval) / n < 1e-2,
                  info = paste(kind, bwmethod, ukertype, okertype, "fval per observation"))
    } else {
      expect_true(abs(bw.profile$fval - bw.dense$fval) < 1e-3,
                  info = paste(kind, bwmethod, ukertype, okertype, "fval"))
    }
  }

  run_case("unordered", "cv.ml", "aitchisonaitken", "wangvanryzin", 2001L)
  run_case("unordered", "cv.ls", "liracine", "wangvanryzin", 2002L)
  run_case("ordered", "cv.ml", "aitchisonaitken", "racineliyan", 2003L)
  run_case("ordered", "cv.ls", "aitchisonaitken", "liracine", 2004L)
  run_case("mixed", "cv.ml", "liracine", "racineliyan", 2005L)
  run_case("mixed", "cv.ls", "liracine", "liracine", 2006L)
})
