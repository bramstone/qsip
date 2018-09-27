#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities have not been calculated
#'   beforehand, \code{calc_d_wad} will compute those first.
#'
#' @return \code{calc_d_wad} adds two S4 Matrix objects to the \code{data@@qsip@@.Data} slot, one for differences
#'   in weighted average density, and the other for weighted average density values of light treatments only (to be used in
#'   future calculations). The row and column specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}},
#'   meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix objects
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
  if(length(data@qsip@iso_trt)==0) stop('Must specify treatment and controls with iso_trt')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # manipulate data matrix and calculate
  # split by replicate groups, but keep track of light and heavy fractions
  iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  # Drop any rows (probably NA) that don't appear in ft rownames
  keep_rows <- iso_group$replicate %in% rownames(ft)
  if(sum(!keep_rows) > 0) {
    warning('Dropping group(s) ',
            as.character(iso_group$replicate[!keep_rows]),
            ' from calculation', call.=FALSE)
  }
  iso_group <- iso_group[keep_rows,]
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # calculate average WAD per taxa for each replicate group
  ft <- base::lapply(ft, colmeans, na.rm=T)
  # create a new list to add results of mean WAD difference into
  d_ft <- as.list(rep(0, nlevels(iso_group$grouping)))
  d_ft <- base::lapply(d_ft, matrix,
                 rep(0, phyloseq::ntaxa(data)),
                 nrow=1,
                 ncol=phyloseq::ntaxa(data))
  names(d_ft) <- levels(iso_group$grouping)
  # For each repliate group: identify which elements of ft are light and which are heavy, then get difference
  # ALTERNATIVE, CALCULATE AVERAGE WAD_LIGHT VALUES EXPERIMENT-WIDE
  iso_group2 <- unique(iso_group[,!names(iso_group) %in% 'replicate']) # only get unique elements to match levels in ft
  for(i in 1:length(d_ft)) {
    which_light <- which(as.numeric(iso_group2$grouping)==i &
                           as.numeric(iso_group2$iso)==1)
    which_heavy <- which(as.numeric(iso_group2$grouping)==i &
                           as.numeric(iso_group2$iso)==2)
    light <- ft[[which_light]]
    heavy <- ft[[which_heavy]]
    d_ft[[i]] <- heavy - light
  }
  # organize and add new data as S4 matrix
  data <- collate_results(data, d_ft, 'd_wad')
  # return weighted average densities of light calcs only
  wl <- ft[[as.numeric(iso_group2$iso)==1]]
  data <- collate_results(data, wl, 'wad_light')
  return(data)
}
