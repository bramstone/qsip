#' Subsetting functions
#' @inheritParams base::Extract
#' @exportMethod

setMethod("[[", "qsip", function(x, i, ...) {
  match_name <- which(attributes(x)$names==i)
  if(length(match_name)==0) return(NULL)
  as.list(x@.Data)[[match_name]]
})
