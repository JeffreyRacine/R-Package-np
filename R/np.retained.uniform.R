# A retained constructor context is not a fresh user request for kernel order.
# Callers must first establish retained-object provenance; public constructors
# and their ordinary fresh-input advisories remain unchanged.
.np_retained_uniform_constructor <- function(constructor, ckertype, ...) {
  withCallingHandlers(constructor(ckertype = ckertype, ...),
    warning = function(w) {
      if (identical(ckertype, "uniform") &&
          identical(conditionMessage(w), unname(.np_io_prefix_text(
            "ignoring kernel order specified with uniform kernel type"))))
        invokeRestart("muffleWarning")
    })
}
