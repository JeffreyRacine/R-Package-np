rank_support_fixture <- function() {
  set.seed(4711)
  x <- data.frame(x = runif(40, .05, .95), z = runif(40, .05, .95))
  list(x = x, y = sin(4*x$x) + x$z^2 + rnorm(40, sd = .2),
       ex = x[c(3, 9, 14, 22, 31), ])
}

test_that("structurally sparse LL hats and fits share the canonical ridge", {
  pkg <- environmentName(environment(npreg))
  f <- rank_support_fixture()
  bw <- suppressWarnings(npregbw(xdat = f$x, ydat = f$y, bws = c(11, 11),
    bandwidth.compute = FALSE, bwtype = "generalized_nn", ckertype = "uniform",
    regtype = "ll"))
  kw <- getFromNamespace(".np_kernel_weights_direct", pkg)(
    bws = bw, txdat = f$x, exdat = f$ex, bandwidth.divide = TRUE)
  basis <- getFromNamespace("W.lp", pkg)
  z <- basis(f$x, degree = c(1L, 1L), basis = "glp", bernstein.basis = FALSE)
  d <- basis(f$x, exdat = f$ex, degree = c(1L, 1L), basis = "glp",
             bernstein.basis = FALSE, gradient.vec = c(1L, 0L))
  expect_identical(sum(kw[, 1L] != 0), 2L)
  a <- crossprod(z, z * kw[, 1L])
  delta <- max(abs(diag(a))) / nrow(z)
  u <- solve(a + diag(delta, ncol(z)), d[1L, ])
  u[1L] <- u[1L] * (1 + delta / a[1L, 1L])
  oracle <- kw[, 1L] * drop(z %*% u)
  old <- getOption("matprod")
  on.exit(options(matprod = old), add = TRUE)
  for (mode in c("default", "blas", "internal")) {
    options(matprod = mode)
    h <- npreghat(bw, txdat = f$x, exdat = f$ex, s = c(1L, 0L),
                  output = "matrix")
    expect_equal(as.double(h[1L, ]), oracle, tolerance = 2e-12, info = mode)
    for (y in list(f$y, cos(f$x$x), rep(1, 40))) {
      g <- gradients(npreg(bw, txdat = f$x, tydat = y, exdat = f$ex,
                            gradients = TRUE))[, 1L]
      expect_equal(drop(h %*% y), g, tolerance = 2e-10, info = mode)
    }
  }
})

test_that("donor certificates count rows rather than bootstrap multiplicity", {
  pkg <- environmentName(environment(npreg))
  count <- getFromNamespace(".np_inid_lp_rank_bounds", pkg)
  project <- getFromNamespace(".np_inid_lp_batch_project", pkg)
  w <- c(0, 1, -1, 0, 2)
  counts <- cbind(c(0, 8, 0, 0, 0), c(0, 1, 1, 0, 6), rep(1, 5))
  expect_identical(count(w, 3L, counts), c(1L, 3L, 3L))
  for (k in 0:4)
    expect_identical(count(c(rep(1, k), rep(0, 5-k)), 3L), as.integer(min(3, k)))
  z <- cbind(1, seq(-1, 1, length.out = 5), seq(-1, 1, length.out = 5)^2)
  y <- sin(seq_len(5))
  cw <- w * counts[, 1L]
  a <- crossprod(z, z*cw)
  rhs <- drop(crossprod(z, y*cw))
  packed <- unlist(lapply(1:3, function(j) a[j, j:3]))
  result <- project(matrix(packed, 1), matrix(rhs, 1), c(1, .2, .04),
                    represented.mass = 8, diagnostics = TRUE, rank.bounds = 1L)
  expect_identical(result$ridge_steps, 1L)
  delta <- max(abs(diag(a))) / 8
  rhs[1L] <- rhs[1L] * (1 + delta/a[1L, 1L])
  expected <- sum(c(1, .2, .04) * solve(a + diag(delta, 3), rhs))
  expect_equal(result$values[1L, 1L], expected, tolerance = 2e-13)
  expect_error(project(matrix(packed, 1), matrix(rhs, 1), c(1, 0, 0),
                       represented.mass = 8, rank.bounds = -1L), "rank_bounds")
})
test_that("higher-degree sparse rows use the same response and adjoint policy", {
  pkg <- environmentName(environment(npreg))
  x <- data.frame(x = c(-.1, 0, .1, seq(1, 3, length.out=9)),
                  z = c(.1, -.1, 0, 2.3, 1.1, 2.7, 1.4, 3, 1.8, 2.1, 1.3, 2.6))
  y <- sin(x$x) + cos(x$z)
  ex <- data.frame(x = 0, z = 0)
  for (basis.name in c("glp", "additive", "tensor"))
  for (bernstein in c(FALSE, TRUE)) {
    bw <- npregbw(xdat=x, ydat=y, bws=c(.3,.3), regtype="lp",
      degree=c(2L,2L), degree.select="manual", basis=basis.name,
      bernstein.basis=bernstein,
      ckertype="uniform", bandwidth.compute=FALSE)
    kw <- getFromNamespace(".np_kernel_weights_direct",pkg)(
      bws=bw, txdat=x, exdat=ex, bandwidth.divide=TRUE)[,1L]
    basis <- getFromNamespace("W.lp",pkg)
    z <- basis(x,degree=c(2L,2L),basis=basis.name,bernstein.basis=bernstein)
    expect_lt(sum(kw!=0),ncol(z))
    a <- crossprod(z,z*kw)
    delta <- max(abs(diag(a)))/nrow(z)
    for (s in list(c(0L,0L),c(1L,0L),c(0L,1L))) {
      d <- basis(x,exdat=ex,degree=c(2L,2L),basis=basis.name,
        bernstein.basis=bernstein,gradient.vec=if(sum(s)) s else NULL)
      u <- solve(a+diag(delta,ncol(z)),d[1L,])
      u[1L] <- u[1L]*(1+delta/a[1L,1L])
      oracle <- kw*drop(z%*%u)
      h <- npreghat(bw,txdat=x,exdat=ex,s=s,output="matrix")
      expect_equal(as.double(h),oracle,tolerance=2e-11)
      fit <- npreg(bw,txdat=x,tydat=y,exdat=ex,gradients=TRUE)
      answer <- if(sum(s)==0L) fitted(fit) else gradients(fit)[,which(s!=0L)]
      expect_equal(drop(h%*%y),as.double(answer),tolerance=2e-11)
    }
  }
})
test_that("duplicate complete designs cannot inflate the structural rank", {
  pkg <- environmentName(environment(npreg))
  identity <- getFromNamespace(".np_inid_lp_design_identity", pkg)
  count <- getFromNamespace(".np_inid_lp_rank_bounds", pkg)
  W <- rbind(c(1, -0, 0), c(1, 0, -0), c(1, 2, 3), c(1, 2, 3),
             c(1, 4, 9), c(1, 5, 10))
  ids <- identity(W)
  expect_identical(ids[, 1L], c(0L, 0L, 2L, 2L, 4L, 5L))
  weights <- c(0, -1, 1, 1, 0, 0)
  counts <- cbind(rep(1L, 6), c(0L, 3L, 0L, 7L, 0L, 0L), rep(0L, 6))
  for (p in 1:5)
    expect_identical(count(weights, p, counts, ids),
                     as.integer(c(min(p, 2), min(p, 2), 0)))
  expect_null(identity(W[c(1, 3, 5, 6), , drop = FALSE]))
  bad <- ids; bad[1, 2] <- 0L
  expect_error(count(rep(0, 6), 3L, identity = bad), "identity chain")
})

test_that("expanded duplicate rows retain the basis-general ridge oracle", {
  pkg <- environmentName(environment(npreg))
  old <- options(np.messages = FALSE, np.tree = getOption("np.tree"))
  on.exit(options(old), add = TRUE)
  base <- data.frame(x=c(-.1,.1,seq(1,2,length.out=8)),
                     z=c(.15,-.15,2.3,1.1,2.7,1.4,3,1.8,2.1,1.3))
  x <- base[c(rep(1:2,each=6),3:10), ]
  y <- sin(x$x) + cos(x$z)
  ex <- data.frame(x=0,z=0)
  basis <- getFromNamespace("W.lp",pkg)
  for (degree in 1:2) for (bernstein in c(FALSE,TRUE))
  for (tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    b <- suppressWarnings(npregbw(xdat=x,ydat=y,bws=c(.3,.3),
      bandwidth.compute=FALSE,ckertype="uniform",regtype="lp",
      degree=rep(degree,2),degree.select="manual",bernstein.basis=bernstein))
    W <- basis(x,degree=rep(degree,2),basis="glp",bernstein.basis=bernstein)
    E <- basis(x,exdat=ex,degree=rep(degree,2),basis="glp",bernstein.basis=bernstein)
    k <- getFromNamespace(".np_kernel_weights_direct",pkg)(
      bws=b,txdat=x,exdat=ex,bandwidth.divide=TRUE)[,1]
    expect_gt(sum(k != 0),ncol(W))
    expect_equal(nrow(unique(W[k != 0,,drop=FALSE])),2)
    A <- crossprod(W,W*k); delta <- max(abs(diag(A)))/nrow(W)
    rhs <- drop(crossprod(W,y*k)); rhs[1] <- rhs[1]*(1+delta/A[1,1])
    oracle <- drop(E %*% solve(A+diag(delta,ncol(W)),rhs))
    fit <- suppressWarnings(npreg(b,exdat=ex))
    hat <- suppressWarnings(npreghat(b,exdat=ex))
    expect_equal(as.numeric(fitted(fit)),oracle,tolerance=1e-10)
    expect_equal(drop(hat %*% y),oracle,tolerance=1e-10)
  }
})

test_that("conditional rows and large RHS blocks share duplicate-design support", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)
  base <- data.frame(x=c(-.1,.1,seq(1,2,length.out=8)),
    z=c(.15,-.15,2.3,1.1,2.7,1.4,3,1.8,2.1,1.3))
  x <- base[c(rep(1:2,each=6),3:10), ]
  y <- sin(x$x) + cos(x$z)
  b <- suppressWarnings(npregbw(x,y,bws=c(.3,.3),regtype="ll",
    ckertype="uniform",bandwidth.compute=FALSE))
  Y <- outer(y,seq_len(19),function(a,b) sin(a*b))
  H <- suppressWarnings(npreghat(b))
  A <- suppressWarnings(npreghat(b,y=Y,output="apply"))
  expect_equal(as.numeric(A),as.numeric(H%*%Y),tolerance=1e-10)
  ex <- data.frame(x=rep(0,3),z=rep(0,3)); ey <- c(.1,.3,.5)
  for(family in c("npcdens","npcdist")) {
    bwfun <- get(paste0(family,"bw"),asNamespace("np"))
    fitfun <- get(family,asNamespace("np"))
    cb <- suppressWarnings(bwfun(xdat=x,ydat=y,bws=c(.4,.3,.3),
      regtype="ll",cxkertype="uniform",bandwidth.compute=FALSE))
    fit <- suppressWarnings(fitfun(cb,exdat=ex,eydat=ey))
    target <- if(family=="npcdens")
      outer(y,ey,function(a,b) dnorm((b-a)/.4)/.4) else
      outer(y,ey,function(a,b) pnorm((b-a)/.4))
    expected <- vapply(seq_along(ey),function(j)
      as.numeric(fitted(suppressWarnings(npreg(b,tydat=target[,j],
        exdat=ex[j,,drop=FALSE])))),numeric(1))
    expect_lt(max(abs(as.numeric(fitted(fit))-expected)),1e-10)
  }
})
