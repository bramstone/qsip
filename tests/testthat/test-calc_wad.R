#test of calc_wad function using testthat 

library(data.table)

test_that("calc_wad produces expected output", {
  data(example_qsip)
  example_qsip[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
  example_qsip <- freq_filter(example_qsip, 
                            min_freq = 3, 
                            filter_target = 'sampleID',
                            tax_id = 'asv_id')
  wads <- calc_wad(example_qsip, tax_id = "asv_id", sample_id = "sampleID", frac_id = "fraction",
                   frac_dens = "Density.g.ml", frac_abund = "avg_16S_g_soil", rel_abund = "rel_abund",
                   grouping_cols = c("treatment", "isotope", "iso_trt", "Phylum"))
  expect_s3_class(wads, "data.table")
  expect_true(nrow(wads) > 0)
  expect_true(all(c("wad", "wvd", "abund") %in% names(wads)))
})