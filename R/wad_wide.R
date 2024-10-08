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
#'  If \code{iso_trt} is not a factor, then it will be coerced to one and it will
#'  return a message telling the user which value was assigned to which isotope labeling level.
#'  Most users adopt a scheme of using the terms "light" and "label" to define their
#'  treatments but please note that the word "label" comes first in the alphabet and
#'  will be assigned by the function to be the unlabeled or light treatment.
#'  If you wish to avoid this (you almost certainly do), you must specify your factor levels beforehand.
#'
#'
#' @seealso \code{\link{calc_wad}}, \code{\link{calc_excess}}
#'
#' @examples
#'  # Load in example data
#'  data(example_qsip)
#'
#'  # relativize sequence abundances (should be done after taxonomic filtering)
#'  example_qsip[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
#'
#'  # ensure that the "light" treatment is the first factor level in the isotope treatment column
#'  levels(example_qsip$iso_trt)
#'
#'  # calculate weighted average densities
#'  wads <- calc_wad(example_qsip,
#'                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
#'                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
#'                   rel_abund = 'rel_abund',
#'                   grouping_cols = c('treatment', 'isotope', 'iso_trt', 'Phylum'))
#'
#'  # transform to wide format
#'  ww <- wad_wide(wads, tax_id = 'asv_id', sample_id = 'sampleID',
#'                 iso_trt = 'iso_trt', isotope = 'isotope')
#'
#' @export

wad_wide <- function(data, tax_id = c(), sample_id = c(), wads = 'wad',
                     iso_trt = c(), isotope = c(), average_light = TRUE) {
  vars <- list(tax_id, sample_id, iso_trt, isotope)
  if(any(sapply(vars, is.null))) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'isotope addition (for each row "13C" or "15N" or "18O")',
                         'isotope treatment (amended or unamended)')[null_vars],
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
    message('Assigned user value "', light_trt, '" as the unamended or light treatment and user value "',
            levels(test_trt)[!levels(test_trt) %in% light_trt],
            '" as the heavy treatment(s).')
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

