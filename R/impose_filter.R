#' Filters taxa from data set
#'
#' Filters taxa from a feature table matrix based on replicate and fraction frequencies
#'
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param code Optional character vector specifying a particular combination of replicate and fraction frequency to test.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#'
#' @details Some text here
#'
#' @return Some text here
#'
#' @seealso \code{\link{create_filters}}
#'
#' @examples
#' # Filter a tax table
#'
#' @export

impose_filter <- function(replicate=0, fraction=0, code=character(), data) {
  # if !is.null(code) parse and use code
  # extract feature table and convert to matrix with taxa as columns
  ft <- as(data@otu_table, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # make presence-absence
  ft <- ceiling(ft / max(ft))
  storage.mode(ft) <- 'integer'
  # split by replicates
  ft <- split_data(data, ft, data@qsip@rep_id)
  # calculate within replicate (i.e., fraction frequency)
  ft <- lapply(ft, colSums, na.rm=T)
  # apply filter, combine
  ft <- lapply(ft, function(x) ifelse(x >= fraction, 1, 0))
  ft <- do.call(rbind, ft)
  # split by replicate group
  iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  # Drop any rows (probably NA) that don't appear in ft rownames, also drop any rows with no taxa
  keep_rows <- (iso_group$replicate %in% rownames(ft) & rowSums(ft) > 0)
  iso_group <- iso_group[keep_rows,]
  ft <- ft[!is.na(iso_group$iso),]
  iso_group <- iso_group[!is.na(iso_group$iso),]
  iso_group <- iso_group[match(rownames(ft), iso_group$replicate),] # match row order to ft
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # apply filter, combine
  ft <- lapply(ft, colSums, na.rm=T)
  ft <- lapply(ft, function(x) ifelse(x >= replicate, 1, 0))
  ft <- do.call(rbind, ft)
  # add in taxa names and calculate presence or absence after filters
  colnames(ft) <- phyloseq::taxa_names(data)
  output <- colSums(ft)
  output <- ifelse(output > 0, 1, 0)
}
