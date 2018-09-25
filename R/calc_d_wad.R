#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities have not been calculated
#'   beforehand, \code{calc_d_wad} will compute those first.
#'
#' @return \code{calc_d_wad} adds an S4 Matrix class (which more efficiently stores sparse matrix data) to the \code{.Data} slot within
#'   the \code{qsip} slot of the object of differences in weighted average density values for each taxon at each group of replicates. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate weighted average density differences
#'
#' @export

calc_d_wad <- function(data) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if WAD values don't exist, calculate those first, this will also handle rep_id validity
  if(is.null(data@qsip[['wad']])) data <- calc_wad(data)
  if(length(data@qsip@rep_group)==0) stop('Must specify replicate groupings with rep_group')
  if(length(data@qsip@iso)==0) stop('Must specify treatment and controls with iso')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # manipulate data matrix and calculate
  # split by replicate groups, but keep track of light and heavy fractions
  iso_group <- iso_grouping(data, data@qsip@iso, data@qsip@rep_group)
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # calculate average WAD per taxa for each replicate group
  ft <- base::apply(ft, 2, mean, na.rm=T)
  # create a new list to add results of mean WAD difference into
  d_ft <- as.list(rep(0,nlevels(iso_grouping$grouping)))
  d_ft <- base::lapply(d_ft, matrix,
                 rep(0, phyloseq::ntaxa(data)),
                 nrow=1,
                 ncol=phyloseq::ntaxa(data))
  names(d_ft) <- levels(iso_group$grouping)
  # For each repliate group: identify which elements of ft are light and which are heavy, then get difference
  for(i in 1:length(d_ft)) {
    which_light <- which(as.numeric(iso_group$grouping)==i &
                           as.numeric(iso_group$iso)==1)
    which_heavy <- which(as.numeric(iso_group$grouping)==i &
                           as.numeric(iso_group$iso)==2)
    light <- ft[[which_light]]
    heavy <- ft[[which_heavy]]
    d_ft[[i]] <- heavy - light
  }
  # organize and add new data as S4 matrix
  data <- collate_results(data, ft, 'd_wad')
  return(data)
}
