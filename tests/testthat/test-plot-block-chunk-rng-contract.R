test_that("block samples and RNG do not depend on chunk boundaries", {
  get <- function(name)getFromNamespace(name,"np")
  for(sim in c("fixed","geom")) {
    indices <- get(".np_block_indices_drawer")(n=7L,B=5L,blocklen=3L,sim=sim)
    set.seed(831)
    anchor <- do.call(cbind,lapply(1:5,function(i)indices(i,i)))
    state <- .Random.seed
    for(chunk in c(2L,5L)) {
      set.seed(831)
      actual <- do.call(cbind,lapply(seq.int(1L,5L,by=chunk),function(i)
        indices(i,min(5L,i+chunk-1L))))
      expect_identical(actual,anchor)
      expect_identical(.Random.seed,state)
    }
    counts <- get(".np_block_counts_drawer")(n=7L,B=5L,blocklen=3L,sim=sim)
    set.seed(831)
    expect_identical(counts(1L,5L),vapply(1:5,function(j)
      as.double(tabulate(anchor[,j],nbins=7L)),numeric(7L)))
    expect_identical(.Random.seed,state)
    expect_error(counts(0L,1L),"invalid block bootstrap chunk bounds")
    expect_error(indices(1L,6L),"invalid block bootstrap chunk bounds")
  }
})

test_that("geometric allocation preserves boot's single-replicate RNG", {
  ts.array <- utils::getFromNamespace("ts.array", "boot")
  make.ends <- utils::getFromNamespace("make.ends", "boot")
  for(endcorr in c(FALSE, TRUE)) {
    draw <- getFromNamespace(".np_block_replicate_drawer", "np")(
      n=7L, n.sim=11L, blocklen=2.5, sim="geom", endcorr=endcorr)
    set.seed(519)
    reference <- replicate(3L, {
      z <- ts.array(n=7L, n.sim=11L, R=1L, l=2.5,
                    sim="geom", endcorr=endcorr)
      ends <- cbind(as.vector(z$starts), as.vector(z$lengths))
      indices <- apply(ends, 1L, make.ends, 7L)
      as.integer(unlist(indices))[seq_len(11L)]
    }, simplify=FALSE)
    state <- .Random.seed
    set.seed(519)
    expect_identical(replicate(3L, draw(), simplify=FALSE), reference)
    expect_identical(.Random.seed, state)
  }
})
