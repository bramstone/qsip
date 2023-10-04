#' Wide-format for WAD values
#'
#' Transforms weighted average density values into wide-format
#'
#' @param data Data as a long-format data.table where each row represents a taxonomic feature within a single fraction.
#'  Typically, this is the output from the \code{calc_wad} function.
#' @param tax_id Column name specifying unique identifier for each taxonomic feature.
#' @param sample_id Column name specifying unique identifier for each replicate.
#' @param wads Column name specifying weighted average density values.
#' @param iso_trt Column name specifying a two-level categorical column indicating whether a sample has been amended with a stable isotope (i.e., is "heavy") or if
#'  isotopic composition is at natural abundance (i.e., "light").
#'  Any terms may be applied but care should be taken for these values.
#'  If supplied as a factor, \code{calc_excess} will take the lowest level as the "light" treatment and the higher
#'  level as the "heavy" treatment.
#'  Alternatively, if supplied as a character, \code{calc_excess} will coerce the column to a factor and with the default behavior wherein the
#'  first value in alphabetical order will be assumed to be the lowest factor level (i.e. the "light" treatment).
#' @param average_light Whether to average the "light" WAD values.
#'
#' @return \code{wad_wide} returns a wide-transformed data.table where each row represents a taxonomic feature within a single replicate.
#'  If averaging unamended ("light") values is not performed, each taxon will have one numeric value and one \code{NA}
#'  value based on whether it was present in an unamended "light" sample or an amended "heavy" sample.
#'
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_excess}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Transform to wide format
#'
#'
#' @export

wad_wide <- function(data, tax_id = c(), sample_id = c(), wads = 'wad',
                     iso_trt = c(), isotope = c(), average_light = TRUE) {
  vars <- list(tax_id, sample_id, iso_trt, isotope, wads)
  if(any(sapply(vars, is.null))) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'isotope addition (for each row "13C" or "15N" or "18O")',
                         'isotope treatment (amended or unamended)', 'WAD values')[null_vars],
                       sep = ',')
    stop("Must supply the following columns: ", null_vars)
  }
  if(any(sapply(vars, function(x) !exists(x, data)))) {
    missing_vars <- which(sapply(vars, function(x) !exists(x, data)))
    missing_vars <- paste(vars[missing_vars], sep = ',')
    stop("Missing the following column(s) in supplied data: ", missing_vars)
  }
  ww <- copy(data)
  setnames(data, old = wads, new = 'wad')
  # make sure iso_trt column is a factor
  if(is.factor(ww[[iso_trt]]) == FALSE || nlevels(ww[[iso_trt]]) > 2) {
    test_trt <- factor(ww[[iso_trt]])
    light_trt <- levels(test_trt)[1]
    message('Assigned ', light_trt, ' as the unamended or "light" treatment and ',
            levels(test_trt)[!levels(test_trt) %in% light_trt],
            ' as the "heavy" treatment(s).')
    ww[[iso_trt]] <- factor(ww[[iso_trt]])
  }
  # rename factor levels
  ww[[iso_trt]] <- factor(ww[[iso_trt]], level = levels(ww[[iso_trt]]), labels = c('light', 'label'))
  # convert to wide format
  all_cols <- setdiff(names(ww), c(iso_trt, wads))
  wide_formula <- paste0(paste(all_cols, collapse = ' + '), ' ~ iso_trt')
  ww <- dcast(ww, as.formula(wide_formula), value.var = 'wad', fill = NA)
  # average light WADs by taxon
  if(average_light) ww[, light := mean(light, na.rm = T), by = tax_id]
  # replace NaN values with NA
  ww[is.nan(light), light := NA]
  return(ww[])
}

