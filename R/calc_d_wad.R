#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
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

calc_d_wad <- function(data, filter=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if WAD values don't exist, calculate those first, this will also handle rep_id validity
  if(is.null(data@qsip[['wad']])) data <- calc_wad(data, filter=filter)
  if(length(data@qsip@rep_group)==0) stop('Must specify replicate groupings with rep_group')
  if(length(data@qsip@iso_trt)==0) stop('Must specify treatment and controls with iso_trt')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  n_taxa <- ncol(ft)
  tax_names <- colnames(ft)
  # manipulate data matrix and calculate
  # split by replicate groups, but keep track of light and heavy fractions
  iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  # Drop any rows (probably NA) that don't appear in ft rownames, also drop any rows with NA for isotope
  keep_rows <- (iso_group$replicate %in% rownames(ft) & !is.na(iso_group$iso))
  if(sum(!keep_rows) > 0) {
    warning('Dropping group(s): ',
            paste(as.character(iso_group$replicate[!keep_rows]), collapse=', '),
            ' - from calculation', call.=FALSE)
  }
  iso_group <- iso_group[iso_group$replicate %in% rownames(ft),]
  ft <- ft[!is.na(iso_group$iso),]
  iso_group <- iso_group[!is.na(iso_group$iso),]
  iso_group <- iso_group[match(rownames(ft), iso_group$replicate),] # match row orders to ft
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # WAD values of 0 indicate no taxa presence in that replicate, convert to NA
  # so that mean WAD values are not pulled down
  ft <- base::lapply(ft, function(x) {x[x==0] <- NA; x})
  # calculate average WAD per taxa for each replicate group
  ft <- base::lapply(ft, colMeans, na.rm=T)
  # remove any NaNs resulting from when a taxon is missing in all replicates
  ft <- base::lapply(ft, function(x) {x[is.nan(x)] <- NA; x})
  # create a new list to add results of mean WAD difference into
  d_ft <- as.list(rep(0, nlevels(iso_group$grouping)))
  d_ft <- base::lapply(d_ft, matrix,
                       rep(0, n_taxa),
                       nrow=1,
                       ncol=n_taxa)
  names(d_ft) <- levels(iso_group$grouping)
  # For each repliate group: identify which elements of ft are light and which are heavy, then get difference
  iso_group2 <- unique(iso_group[,!names(iso_group) %in% 'replicate']) # only get unique elements to match levels in ft
  for(i in 1:length(d_ft)) {
    # use numbers to reference non-labeled additions since they're element agnostic
    # any NA values result here when a taxa is completely missing from a heavy or light treatment in a replicate group
    which_light <- which(as.numeric(iso_group2$grouping)==i &
                           as.numeric(iso_group2$iso)==1)
    which_heavy <- which(as.numeric(iso_group2$grouping)==i &
                           as.numeric(iso_group2$iso)==2)
    light <- ft[[which_light]]
    heavy <- ft[[which_heavy]]
    d_ft[[i]] <- heavy - light
  }
  # organize and add new data as S4 matrix
  data <- collate_results(data, d_ft, 'd_wad', tax_names=tax_names, sparse=TRUE)
  # return weighted average densities of light calcs only
  iso_group2 <- iso_group2[match(names(ft), iso_group2$interaction),] # match row order to ft
  wl <- ft[which(as.numeric(iso_group2$iso)==1)]
  data <- collate_results(data, wl, 'wad_light', tax_names=tax_names, sparse=TRUE)
  return(data)
}
