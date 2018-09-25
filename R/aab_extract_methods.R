#' Subsetting functions
#' @exportMethod

setMethod("[[", "qsip", function(x, i, ...) {
  list(x)[[which(attributes(x)$names==i)]]
})
