#test of calc_pop function using testthat 

library(data.table)

test_that("calc_pop produces expected non-bootstrapped output", {
  # Load and prepare test data
  data(example_qsip)
  
  # Filter out non-target lineages as shown in README
  tx <- unique(example_qsip[, c('asv_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
  unassign <- tx[Kingdom == 'Unassigned', asv_id]
  euk <- tx[Kingdom == 'Eukaryota', asv_id]
  arch <- tx[Kingdom == 'Archaea', asv_id]
  mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')
                    ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), asv_id]
  
  example_qsip_filtered <- example_qsip[!asv_id %in% c(arch, mito_chloro, euk, unassign)]
  
  # Relativize sequence abundances
  example_qsip_filtered[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
  
  # Calculate WADs first (required input for calc_pop)
  wads <- calc_wad(example_qsip_filtered,
                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
                   rel_abund = 'rel_abund',
                   grouping_cols = c('treatment', 'isotope', 'iso_trt', 'Phylum', 'timepoint'))

  # create time-0 data, ensure abundance column is present with same name as WADs data
  t0_samples <- example_qsip[timepoint == 0]
  t0_samples[, abund := avg_16S_g_soil * (seq_abund / sum(seq_abund)), by = sampleID]
  
  # Test calc_pop non-bootstrapped
  rates <- calc_pop(wads, t0_samples, tax_id = 'asv_id', sample_id = 'sampleID',
                    iso_trt = 'iso_trt', isotope = 'isotope', time = 'timepoint')

# Tests..................................................
  
  # Basic structure tests
  expect_s3_class(rates, "data.table")
  expect_true(nrow(rates) > 0)
  
  # Check for expected columns
  expected_cols <- c("growth", "mortality")
  expect_true(all(expected_cols %in% names(rates)))
  
  # Check that rate values are numeric
  expect_true(is.numeric(rates$growth))
  expect_true(is.numeric(rates$mortality))
  
  # Check that we have the expected grouping columns
  expected_grouping_cols <- c('asv_id', 'treatment', 'Phylum')
  expect_true(all(expected_grouping_cols %in% names(rates)))
  
  # Check that rates are not all NA (should have some valid calculations)
  expect_true(sum(!is.na(rates$growth)) > 0)
  expect_true(sum(!is.na(rates$mortality)) > 0)
  
  # Check that rates are reasonable (not extremely large values)
  finite_growth <- rates$growth[is.finite(rates$growth)]
  finite_mortality <- rates$mortality[is.finite(rates$mortality)]
  if(length(finite_growth) > 0) {
    expect_true(all(abs(finite_growth) < 1000))  # Reasonable range check
  }
  if(length(finite_mortality) > 0) {
    expect_true(all(abs(finite_mortality) < 1000))  # Reasonable range check
  }
})

test_that("calc_pop produces expected bootstrapped output", {
  # Load and prepare test data
  data(example_qsip)
  
  # Filter out non-target lineages
  tx <- unique(example_qsip[, c('asv_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
  unassign <- tx[Kingdom == 'Unassigned', asv_id]
  euk <- tx[Kingdom == 'Eukaryota', asv_id]
  arch <- tx[Kingdom == 'Archaea', asv_id]
  mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')
                    ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), asv_id]
  
  example_qsip_filtered <- example_qsip[!asv_id %in% c(arch, mito_chloro, euk, unassign)]
  
  # Relativize sequence abundances
  example_qsip_filtered[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
  
  # Calculate WADs first
  wads <- calc_wad(example_qsip_filtered,
                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
                   rel_abund = 'rel_abund',
                   grouping_cols = c('treatment', 'isotope', 'iso_trt', 'Phylum', 'timepoint'))

  # create time-0 data, ensure abundance column is present with same name as WADs data
  t0_samples <- example_qsip[timepoint == 0]
  t0_samples[, abund := avg_16S_g_soil * (seq_abund / sum(seq_abund)), by = sampleID]
  
  # Test calc_pop with bootstrapping (using fewer iterations for speed)
  rates_boot <- calc_pop(wads, t0_samples, tax_id = 'asv_id', sample_id = 'sampleID',
                         iso_trt = 'iso_trt', isotope = 'isotope', time = 'timepoint',
                         bootstrap = TRUE, iters = 99,
                         grouping_cols = c('treatment', 'Phylum'))
  
# Tests..................................................
  
  # Basic structure tests
  expect_s3_class(rates_boot, "data.table")
  expect_true(nrow(rates_boot) > 0)
  
  # Check for bootstrap-specific columns
  bootstrap_cols <- c("2.5%", "50%", "97.5%", "p_val")
  expect_true(all(bootstrap_cols %in% names(rates_boot)))
  
  # Check for rate column that differentiates growth from mortality
  expect_true("rate" %in% names(rates_boot))
  expect_true(all(c("growth", "mortality") %in% unique(rates_boot$rate)))
  
  # Check that bootstrap columns are numeric
  expect_true(is.numeric(rates_boot$`2.5%`))
  expect_true(is.numeric(rates_boot$`50%`))
  expect_true(is.numeric(rates_boot$`97.5%`))
  expect_true(is.numeric(rates_boot$p_val))
  
  # Check that confidence intervals make sense (2.5% <= 50% <= 97.5%) where not NA
  valid_rows <- !is.na(rates_boot$`2.5%`) & !is.na(rates_boot$`50%`) & !is.na(rates_boot$`97.5%`)
  if(sum(valid_rows) > 0) {
    expect_true(all(rates_boot$`2.5%`[valid_rows] <= rates_boot$`50%`[valid_rows]))
    expect_true(all(rates_boot$`50%`[valid_rows] <= rates_boot$`97.5%`[valid_rows]))
  }
  
  # Check that p_val is between 0 and 1
  expect_true(all(rates_boot$p_val >= 0 & rates_boot$p_val <= 1, na.rm = TRUE))
  
  # Check that we have expected grouping columns
  expected_cols <- c('asv_id', 'treatment', 'Phylum')
  expect_true(all(expected_cols %in% names(rates_boot)))
})

test_that("calc_pop handles different mu values", {
  # Load and prepare test data
  data(example_qsip)
  
  # Filter out non-target lineages
  tx <- unique(example_qsip[, c('asv_id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
  unassign <- tx[Kingdom == 'Unassigned', asv_id]
  euk <- tx[Kingdom == 'Eukaryota', asv_id]
  arch <- tx[Kingdom == 'Archaea', asv_id]
  mito_chloro <- tx[, tax_string := paste(Phylum, Class, Order, Family, Genus, sep = ';')
                    ][grepl('mitochond|chloroplast', tax_string, ignore.case = T), asv_id]
  
  example_qsip_filtered <- example_qsip[!asv_id %in% c(arch, mito_chloro, euk, unassign)]
  
  # Relativize sequence abundances
  example_qsip_filtered[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
  
  # Calculate WADs first
  wads <- calc_wad(example_qsip_filtered,
                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
                   rel_abund = 'rel_abund',
                   grouping_cols = c('treatment', 'isotope', 'iso_trt', 'Phylum', 'timepoint'))

  # create time-0 data, ensure abundance column is present with same name as WADs data
  t0_samples <- example_qsip[timepoint == 0]
  t0_samples[, abund := avg_16S_g_soil * (seq_abund / sum(seq_abund)), by = sampleID]
  
  # Test with default mu (0.6)
  rates_default <- calc_pop(wads, t0_samples, tax_id = 'asv_id', sample_id = 'sampleID',
                            iso_trt = 'iso_trt', isotope = 'isotope', time = 'timepoint')
  
  # Test with different mu value (0.8)
  rates_mu08 <- calc_pop(wads, t0_samples, tax_id = 'asv_id', sample_id = 'sampleID',
                         iso_trt = 'iso_trt', isotope = 'isotope', time = 'timepoint', mu = 0.8)

  # Tests..................................................
  
  # Both should produce valid outputs
  expect_s3_class(rates_default, "data.table")
  expect_s3_class(rates_mu08, "data.table")
  
  # Results should be different (unless all rates are NA)
  if(sum(!is.na(rates_default$growth)) > 0 && sum(!is.na(rates_mu08$growth)) > 0) {
    expect_false(identical(rates_default$growth, rates_mu08$growth))
  }
})

test_that("calc_pop throws appropriate errors for missing parameters", {
  data(example_qsip)
  
  # Prepare minimal dataset
  example_qsip[, rel_abund := seq_abund / sum(seq_abund), by = sampleID]
  wads <- calc_wad(example_qsip,
                   tax_id = 'asv_id', sample_id = 'sampleID', frac_id = 'fraction',
                   frac_dens = 'Density.g.ml', frac_abund = 'avg_16S_g_soil',
                   rel_abund = 'rel_abund',
                   grouping_cols = c('timepoint'))

  # Tests..................................................
  
  # Test missing tax_id
  expect_error(calc_pop(wads, sample_id = 'sampleID', iso_trt = 'iso_trt', time = 'timepoint'))
  
  # Test missing sample_id
  expect_error(calc_pop(wads, tax_id = 'asv_id', iso_trt = 'iso_trt', time = 'timepoint'))
  
  # Test missing iso_trt
  expect_error(calc_pop(wads, tax_id = 'asv_id', sample_id = 'sampleID', time = 'timepoint'))
  
  # Test missing time
  expect_error(calc_pop(wads, tax_id = 'asv_id', sample_id = 'sampleID', iso_trt = 'iso_trt'))
})