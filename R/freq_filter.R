#' Apply frequency filter
#'
#' Filter taxonomic features in data based on frequency
#'
#' @param data data.table object in long-format where each row represents a sequence feature in a given sequenced sample.
#' @param min_freq Numeric value specifying the minimum acceptable frequency for a taxonomic feature within all replicates.
#' @param filter_target Column name or group of columns specifying the ID variable for each sample-fraction.
#'  Column(s) should be able to uniquely differentiate every fraction from every replicate.
#' @param tax_id Column name specifying a unique identifier for each sequence feature.
#'
#' @details There is no strong agreement among qSIP-users on the proper threshold to set for removing rare features.
#'   Rare and infrequent taxa produce noise in the data, making it hard to discern quality.
#'   The one guiding principle that there may be agreement on is that it’s best to set minimum filters at first -- to be as inclusive as possible –- 
#'   and intensify filters as needed to reduce noise.
#'
#'   One or more columns may be specified in the filter_target parameter, allowing for frequency filtering across treatment groups.
#'   My recommendation in this case is to make sure that all non-labeled samples (those where isotopic composition is at natural abundance)
#'   grouped together so that non-labeled buoyant density estimates may be made with as many occurrences as possible for each taxon.
#'
#'
#' @return Returns filtered data table with taxa above specified frequency thresholds
#'
#' @seealso \code{\link{seq_summary}}
#'
#' @examples
#' data(example_data)
#'
#' # initial sequence and ASV count? 
#' seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' # Remove taxa that occur in fewer than 3 fractions in any given replicate
#' dat <- freq_filter(dat, 3, 'sample_id', 'taxon_id')
#'
#' seq_summary(dat, 'seq_abund', 'taxon_id')
#'
#' @export

freq_filter <- function(data, min_freq = 0, filter_target = c(), tax_id = c()) {
  if(is.null(filter_target) || is.null(tax_id)) stop("Must supply columns for sample ID and taxon ID")
  if(length(tax_id) > 1) stop("Must supply single column name")
  if(length(min_freq) > 1) stop("Must supply single frequency threshold value")
  return(data[, freq := .N, by = c(filter_target, tax_id)
              ][freq >= min_freq
                ][, freq := NULL])
}
