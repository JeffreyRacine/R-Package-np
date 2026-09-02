for (nn.type in c("generalized_nn", "adaptive_nn")) {
  for (nn.family in c("regression", "density", "distribution",
                      "conditional-density", "conditional-distribution")) {
    test_that(paste(nn.family, nn.type, "search returns a raw-valid endpoint"), {
      old <- options(np.messages = FALSE, np.tree = FALSE)
      on.exit(options(old), add = TRUE)
      n <- 24L
      continuous <- data.frame(x = c(rep(0, 16L), seq_len(8L)))
      mixed <- data.frame(u = factor(rep(c("a", "b"), n / 2L)), continuous)
      y <- 1 + seq_len(n) / n + sin(seq_len(n)) / 10
      controls <- list(bwtype = nn.type, nmulti = 1L,
                        powell.remin = FALSE, bwsolver = "powell")
      if (nn.family == "regression") {
        bw <- do.call(npregbw, c(list(xdat = mixed, ydat = y, regtype = "lc"), controls))
        raw <- np:::.npregbw_eval_only(mixed, y, bw,
                                       invalid.penalty = "dbmax")$objective
      } else if (nn.family == "density") {
        bw <- do.call(npudensbw, c(list(dat = mixed), controls))
        raw <- np:::npudensbw.bandwidth(
          dat = mixed, bws = bw, eval.only = TRUE, bandwidth.compute = TRUE,
          nmulti = 1L, invalid.penalty = "dbmax"
        )$fval
      } else if (nn.family == "distribution") {
        bw <- do.call(npudistbw, c(list(dat = continuous), controls))
        raw <- np:::npudistbw.dbandwidth(
          dat = continuous, bws = bw, eval.only = TRUE, bandwidth.compute = TRUE,
          nmulti = 1L, invalid.penalty = "dbmax"
        )$fval
      } else if (nn.family == "conditional-density") {
        bw <- do.call(npcdensbw, c(list(xdat = mixed, ydat = data.frame(y)), controls))
        raw <- np:::.npcdensbw_eval_only(mixed, data.frame(y), bw,
                                         invalid.penalty = "dbmax")$objective
      } else {
        bw <- do.call(npcdistbw, c(list(xdat = mixed, ydat = data.frame(y)), controls))
        raw <- np:::.npcdistbw_eval_only(mixed, data.frame(y), bws = bw,
                                         invalid.penalty = "dbmax")$objective
      }
      expect_length(raw, 1L)
      expect_true(is.finite(raw) && abs(raw) < .Machine$double.xmax)
      expect_identical(as.numeric(raw), as.numeric(bw$fval))
    })
  }
}
