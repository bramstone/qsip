#' Matrix to long-form data
#'
#' Converts sample-by-feature matrix to long-form data.table
#'
#' @param data An matrix of abundances, relativized or not, of taxonomic features across samples.
#'  \code{data} may be an S3 or S4 class from the \code{Matrix} package but must be coercible to a \code{TsparseMatrix}.
#'
#' @details This function is included as a convenience for quickly transforming wide-format feature
#'  tables into long-format.
#'
#'  Sometimes operations involving the \code{Matrix}-package classes begin to fail without explanation.
#'  When this happens, it is prudent to restart R in a fresh session.
#'
#' @return Returns a long-form data.table object where each row represents the abundance of a taxonomic feature
#'  in a sample.
#'
#' @examples to go here
#'
#' @export


mtl <- function(data) {
  if(is.null(colnames(data))) colnames(data) <- paste('col', 1:ncol(data))
  if(is.null(rownames(data))) rownames(data) <- paste('row', 1:nrow(data))
  y <- as(x, 'TsparseMatrix')
  return(data.table(rownames = y@Dimnames[[1]][y@i + 1],
                    colnames = y@Dimnames[[2]][y@j + 1],
                    value = y@x))
}
