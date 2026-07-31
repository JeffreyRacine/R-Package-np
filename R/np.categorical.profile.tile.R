.np_categorical_profile_kernel_tile <- function(train.unordered,
                                                train.ordered,
                                                eval.unordered,
                                                eval.ordered,
                                                kernel.unordered,
                                                kernel.ordered,
                                                operator,
                                                lambda,
                                                ncat,
                                                supports,
                                                eval.start = 1L,
                                                eval.count = nrow(eval.unordered)) {
  .Call(
    C_np_categorical_profile_kernel_tile,
    train.unordered,
    train.ordered,
    eval.unordered,
    eval.ordered,
    as.integer(kernel.unordered),
    as.integer(kernel.ordered),
    as.integer(operator),
    as.double(lambda),
    as.integer(ncat),
    supports,
    as.integer(eval.start),
    as.integer(eval.count)
  )
}
