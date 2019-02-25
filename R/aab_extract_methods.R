#' Subsetting functions
#' @inheritParams base::Extract
#' @exportMethod

setMethod("[[", "qsip", function(x, i, ...) {
  match_name <- which(attributes(x)$names==i)
  if(length(match_name)==0) return(NULL)
  y <- as.list(x@.Data)[[match_name]]
  if(!is.null(dim(y))) {
    y <- (as(y, 'matrix'))
  }
  y[y==0] <- NA
  return(y)
})
