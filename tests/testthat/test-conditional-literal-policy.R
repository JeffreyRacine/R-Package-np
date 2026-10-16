# These tests protect the per-namespace specialization and a literal kernel
# oracle. MPI endpoint ownership and installed/source parity are separate gates.
conditional_literal_ns <- function() asNamespace("npRmpi")

test_that("literal native substitution preserves NULL and missing AST nodes", {
  rewrite <- get(".np_conditional_literal_rewrite",conditional_literal_ns())
  original <- quote(function(x=NULL, missing.argument, ...) {
    value <- x[, , drop=FALSE]
    list(NULL, value, .Call("C_np_kernelsum", NULL, x, PACKAGE="np"))
  })
  expected <- quote(function(x=NULL, missing.argument, ...) {
    value <- x[, , drop=FALSE]
    list(NULL, value, .Call(.np_literal_C_np_kernelsum, NULL, x, PACKAGE="np"))
  })
  got <- rewrite(original,"C_np_kernelsum")
  expect_identical(got$count,1L)
  expect_identical(got$body,expected)
  expect_identical(rewrite(quote(list(NULL,x[,],function(a=NULL,b)NULL)),"absent")$body,
                   quote(list(NULL,x[,],function(a=NULL,b)NULL)))
  expect_identical(rewrite(quote(.Call("other",NULL)),"C_np_kernelsum")$count,0L)
  expect_identical(rewrite(quote(list(.Call("C_np_kernelsum"),.Call("C_np_kernelsum"))),
                           "C_np_kernelsum")$count,2L)
})

test_that("conditional literal closures preserve the original generic closures", {
  ns <- conditional_literal_ns()
  holder <- get(".np_conditional_literal_holder", ns)
  lex <- holder$lex
  expect_true(environmentIsLocked(holder))
  expect_true(environmentIsLocked(lex))
  expect_true(bindingIsLocked("lex", holder))
  targets <- c(
    ".np_regression_direct", ".npreghat_exact_matrix_from_core",
    ".npreghat_exact_lp_apply_from_regression_core",
    ".npreghat_exact_lp_matrix_from_regression_core",
    ".npreghat_exact_lp_matrix_from_kernel_weights", ".np_kernel_weights_direct",
    ".npreghat_exact_lc_matrix_from_kernel_weights",
    ".npreghat_exact_lc_derivative_matrix_from_npksum_chunked",
    ".np_local_operator_ksum", ".np_direct_operator_matrix",
    ".np_direct_operator_apply", ".np_exact_operator_matrix",
    ".np_exact_operator_apply", ".npcdhat_make_kernel_matrix",
    ".npcdhat_make_xhat_matrix", "npksum.default",
    ".np_plot_density_derivative_matrix")
  leaves <- c(.np_regression_direct="C_np_regression",
    .npreghat_exact_lp_apply_from_regression_core="C_np_regression_lp_apply_conditional_ctx",
    .np_kernel_weights_direct="C_np_kernelsum",
    .np_local_operator_ksum="C_np_kernelsum",npksum.default="C_np_kernelsum")
  for(name in targets) {
    original <- get(name, ns, inherits=FALSE)
    cloned <- if(name=="npksum.default") lex$.np_conditional_npksum_engine else get(name,lex,inherits=FALSE)
    expect_identical(environment(original),ns,info=name)
    expect_identical(environment(cloned),lex,info=name)
    expect_identical(formals(cloned),formals(original),info=name)
    expect_true(bindingIsLocked(name,lex),info=name)
    if(!name %in% names(leaves)) {
      expect_identical(body(cloned),body(original),info=name)
    } else {
      old <- leaves[[name]]
      count <- 0L
      restore <- function(expr) {
        if(!is.call(expr) && !is.pairlist(expr))return(expr)
        if(is.call(expr) && length(expr)>=2L && is.symbol(expr[[1L]]) &&
           as.character(expr[[1L]]) %in% c(".Call","eval.native") &&
           identical(expr[[2L]],as.name(paste0(".np_literal_",old)))) {
          expr[[2L]] <- old
          count <<- count+1L
        }
        for(i in seq_along(expr)) {
          if(is.null(expr[[i]]) || identical(expr[[i]],quote(expr=)))next
          if(is.call(expr[[i]]) || is.pairlist(expr[[i]]))expr[[i]] <- restore(expr[[i]])
        }
        expr
      }
      expect_identical(restore(body(cloned)),body(original),info=name)
      expect_identical(count,1L,info=name)
      expect_s3_class(get(paste0(".np_literal_",old),lex),"NativeSymbolInfo")
    }
  }
  wrapper <- get("npksum_conditional_literal",ns)
  expect_identical(formals(wrapper),formals(get("npksum.default",ns)))
  expect_identical(lex$npksum.default,wrapper)
  if(identical(unname(getNamespaceName(ns)),"npRmpi")) {
    expect_error(lex$.npRmpi_autodispatch_call(quote(identity(1)),owner.name="npreg"),
                 "unexpected conditional kernel-sum dispatch owner")
  }
})

test_that("literal conditional categorical factors agree with a direct probability oracle", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.largelambda=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=factor(c("a","b","c","a","b","c","a","a","b"),levels=c("a","b","c")))
  y <- data.frame(y=factor(c("u","v","w","u","u","v","w","u","v"),levels=c("u","v","w")))
  ex <- data.frame(x=factor(c("a","b","c"),levels=levels(x$x)))
  ey <- data.frame(y=factor(c("u","v","w"),levels=levels(y$y)))
  lambda <- .664
  bw <- npcdensbw(xdat=x,ydat=y,bws=c(lambda,lambda),
    bandwidth.compute=FALSE,regtype="lc",uxkertype="aitchisonaitken",uykertype="aitchisonaitken")
  literal <- function(tx,ex) ifelse(outer(as.character(ex),as.character(tx),"=="),1-lambda,lambda/2)
  K <- literal(x$x,ex$x);Y <- literal(y$y,ey$y)
  oracle <- rowSums(K*Y)/rowSums(K)
  one <- npcdens(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,se=FALSE)
  uncertainty <- npcdens(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,se=TRUE)
  options(np.largelambda=FALSE)
  two <- npcdens(bws=bw,txdat=x,tydat=y,exdat=ex,eydat=ey,se=FALSE)
  expect_equal(one$condens,oracle,tolerance=2e-12)
  expect_identical(one$condens,two$condens)
  expect_identical(one$condens,uncertainty$condens)
  expect_identical(one[["conditional.kernel.policy",exact=TRUE]],"literal-xy-v1")
  expect_false(getOption("np.largelambda"))
})

test_that("historical conmode values remain readable but mixed uncertainty requires refitting", {
  ns <- conditional_literal_ns()
  draw <- get(".np_plot_conmode_data",ns)
  stamp <- get(".np_conditional_policy_stamp",ns)
  require.policy <- get(".np_conditional_require_replay_policy",ns)
  old <- structure(list(probability.levels=c("a","b"),
    probabilities=cbind(c(.7,.2),c(.3,.8)),probability.errors=matrix(.1,2L,2L),
    xeval=data.frame(x=c(0,1)),xndim=1L,xnames="x",se=TRUE,
    conmode=factor(c("a","b"))),class="conmode")
  retained <- old
  expect_identical(draw(old,errors="none")[[1L]]$probability,c(.7,.2))
  expect_true(all(is.finite(draw(old,errors="asymptotic",band="pointwise")[[1L]]$lower)))
  expect_error(draw(old,errors="asymptotic",band="simultaneous"),"refit the object")
  expect_error(draw(old,errors="bootstrap"),"refit the object")
  expect_error(draw(old,gradients=TRUE,errors="asymptotic"),"refit the object")
  expect_identical(old,retained)
  current <- stamp(old)
  expect_invisible(require.policy(current,"test"))
  expect_identical(current$probabilities,old$probabilities)
  expect_identical(current$probability.errors,old$probability.errors)
  expect_null(old[["conditional.kernel.policy",exact=TRUE]])
})
