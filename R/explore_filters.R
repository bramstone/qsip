#' Explore different filtering levels
#'
#' Generates an S4 \code{filter} class object which stores filtering statistics
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param filters Optional data frame specifying filtering levels to explore.
#'   If blank, \code{explore_filters} will generate statistics for all possible frequencies at both the replicate and fraction level.
#'
#' @details Some text here
#'
#' @return Some text here
#'
#' @seealso \code{\link{create_filters}}, \code{\link{impoose_filter}}
#'
#' @examples
#' # Filter a tax table
#'
#' @export

explore_filters <- function(data, filters) {
  # if no filters provided, will use all
  if(is.null(data@qsip@filter) && is.null(filters)) {
    # number of fractions in data?
    fractions <- table(data@sam_data[[data@qsip@rep_id]])
    # group by isotope treatment, replicate ID, replicate group
    group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    group <- group[match(names(fractions), iso_group$replicate),]
    replicates <- split(fractions, group$interaction) # will drop NA isotope values
    # number of replicates in data?
    relicates <- sapply(replicates, length)
    filters <- create_filters(1:max(replicates), 1:max(fractions))
  }
  # create blank matrix (taxa as rows, number of filter combos as columns)
  output <- matrix(0L, nrow=phyloseq::ntaxa(data),
         ncol=nrow(filters),
         dimnames=list(phyloseq::tax_names(data0),
                       interaction(filters, sep=':')))
  # for every combination, run impose filters, stick output to column i of matrix
  for(i in 1:nrow(filters)) {
    output[,i] <- impose_filter(data,
                       replicate=filters$rep_freq[i],
                       fraction=filters$frac_freq[i])
  }
  # convert to S4 Matrix class
  output <- Matrix::Matrix(output, sparse=T)
  # add other statistics into output
  # NEED TO CREATE NEW S4 CLASS FOR FILTERING STATISTICS
  # NEED TO UPDATE SET_CLASS AND SET_INIT_METHODS SCRIPTS
}
