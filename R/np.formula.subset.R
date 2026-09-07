# Resolve formula subsets before any transport or model-frame replay. Return
# values only: the caller environment is used transiently, never retained.
.np_formula_subset_inputs <- function(data, subset.expr, caller) {
  list(data = data, subset = eval(subset.expr, envir = data, enclos = caller))
}

.np_bws_formula_model_frame <- function(bws, mf.args, data.override = FALSE) {
  call.env <- environment(bws$call)
  if (!data.override && is.environment(call.env) &&
      "data" %in% names(mf.args) && is.language(mf.args[["data"]])) {
    # Resolve only the saved data expression in its existing owner, when
    # model.frame forces it. Formula variables and subset retain their lexical
    # data-mask semantics; neither the stored call nor its environment changes.
    mf.args[["data"]] <- substitute(base::eval(quote(EXPR), envir = OWNER),
                                    list(EXPR = mf.args[["data"]],
                                         OWNER = call.env))
  }
  do.call(stats::model.frame, mf.args,
          envir = environment(mf.args[["formula"]]))
}
