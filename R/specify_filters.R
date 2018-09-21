#' Specify multiple taxa filting options for qSIP calculations
#'
#' Creates a data frame of per-replicate and per-fraction minimum incidence frequencies to retain microbial taxa in the dataset
#'
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param rm_combn Optional character vector specifying any particular combinations of replicate and fraction frequency to remove.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#'
#' @details \code{specify_filter} should be used to specify the desired frequencies of microbial taxa \emph{prior} to calculation of
#'   atom excess fraction / atom percent excess. Filters must also be specified so that diagnostic functions can be run as well.
#'   Because \code{specify_filter} uses \code{\link{base::expand.grid}}, all possible combinations of the first two arguments will
#'   be used to construct the resulting data.frame. If certain combinations are not desirable, the \code{rm_combn} argument can be
#'   supplied to remove them.
#'
#' @return \code{specify_filter} produces a data frame of replicate frequencies and within-replicate-fraction frequencies to be
#'   investigated in downstream analyses.
#'
#' @examples
#' # Only filter on samples
#'  specify_filters(2:3)
#'
#' # Filter on samples and fractions
#' specify_filters(2:3, 10:12)
#'
#' # remove certain combinations
#' specify_filters(2:3, 10:13, c('3:13', '2:10'))
#'
#' @export

specify_filters <- function(replicate=0, fraction=0, rm_combn=character()) {
  filters <- expand.grid(replicate, fraction)
  names(filters) <- c('rep_freq', 'frac_freq')
  filters <- filters[order(filters$rep_freq),]
  if(!missing(rm_combn)) {
    rm_combn <- strsplit(rm_combn, ':')
    rm_combn <- do.call(rbind, rm_combn)
    storage.mode(rm_combn) <- 'integer'
    for(i in 1:nrow(rm_combn)) {
      filters <- filters[filters$rep_freq!=rm_combn[i,1] |
                           filters$frac_freq!=rm_combn[i,2],]
    }
  }
  rownames(filters) <- NULL
  return(filters)
}
