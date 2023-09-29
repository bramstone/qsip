# documentation here

calc_excess <- function(data, tax_id = c(), sample_id = c(), wads = 'wad', 
                        iso_trt = c(), isotope = c(),
                        correction = TRUE, rm_outliers = TRUE, non_grower_prop = 0.1,
                       nat_abund_13C = 0.01111233, nat_abund_15N = 0.003663004, nat_abund_18O = 0.002011429) {
  vars <- list(tax_id, sample_id, iso_trt, isotope, wads))
  if(any(sapply(vars, is.null)) {
    null_vars <- which(sapply(vars, is.null))
    null_vars <- paste(c('taxon IDs', 'sample IDs', 'isotope addition (for each row "13C" or "15N" or "18O")',
                         'isotope treatment (amended or unamended)', 'WAD values')[null_vars],
                       sep = ',')
    stop("Must supply the following columns:", null_vars)
  }
  # re-express the iso_trt column to be either "label" or "light"
  if(!is.factor(iso_trt)) message('Assigned', levels(iso_trt)[1], 'as the unamended or "light" treatment'))
  iso_trt <- as.factor(data$iso_trt)
  iso_trt <- factor(iso_trt, labels = c('light', 'label'))
  # convert to wide format
  all_cols <- setdiff(names(data), c(iso_trt, wads))
  wide_formula <- paste0(paste(all_cols, collapse = ' + '), ' ~ iso_trt')
  data <- dcast(dat, as.formula(wide_formula), value.var = 'wad', fill = NA)
  # average light WADs by taxon
  data[, light := mean(light, na.rm = T), by = tax_id]
