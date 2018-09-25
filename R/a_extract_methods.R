#' Subsetting functions
#' @include set_class.R
#' @exportMethod

setMethod("[[", "qsip", function(x, i, ...) {
  list(x)[[which(attributes(x)$names==i)]]
})
