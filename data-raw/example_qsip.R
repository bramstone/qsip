## code to prepare `example_qsip` dataset for qsip package
# data taken from the following published dataset:
#------------
#   Stone et al. 2022. Nutrients strengthen density dependence of per-capita growth
#   and mortality rates in the soil bacterial community. Oecologia
#   https://doi.org/10.1007/s00442-023-05322-z
# -----------

# data.table package extensively utilized to generate data
library(data.table)

# load data, available on github:
download.file('https://github.com/bramstone/density-dependence-qSIP/raw/main/data/feature_table.tar.gz',
              destfile = '~/Downloads/qsip_example_asv_data.tar.gz')
download.file('https://github.com/bramstone/density-dependence-qSIP/raw/main/data/feature_taxonomy.tar.gz',
              destfile = '~/Downloads/qsip_example_tax_data.tar.gz')
system('mkdir ~/Downloads/qsip_example_data/')
system('tar -xzf ~/Downloads/qsip_example_asv_data.tar.gz -C ~/Downloads/qsip_example_data')
system('tar -xzf ~/Downloads/qsip_example_tax_data.tar.gz -C ~/Downloads/qsip_example_data')

exp <- fread('https://raw.githubusercontent.com/bramstone/density-dependence-qSIP/main/data/experimental_data.csv')
asv <- fread('~/Downloads/qsip_example_data/feature_table.csv')
tax <- fread('~/Downloads/qsip_example_data/taxonomy_table.csv')

# combine
dat <- merge(asv, exp[,c('SampleID', 'isotope', 'treatment',
                         'ecosystem', 'sampleID', 'Density.g.ml',
                         'avg_16S_g_soil', 'timepoint')],
             by = 'SampleID',
             all.y = T)

dat <- merge(dat, tax, by = 'taxonID', all.x = T)

# convert all isotope treatments to either 18O or 16O, remove sample W1_PP_13, it shouldn't be here
dat <- dat[, iso_trt := fifelse(isotope == '18O', 'label', 'light')
           ][timepoint == 0, `:=` (isotope = NA, treatment = NA)
             ][sampleID != 'W1_PP_13']

# identify replicates for abundance matching in the growth/death rate calculations
reps <- unique(dat[, .(sampleID, isotope, treatment, ecosystem, timepoint)
                   ])[order(sampleID)
                      ][, rep := 1:.N, by = .(isotope, treatment, ecosystem, timepoint)]

dat <- merge(dat, reps[, .(sampleID, rep)], all.x = T, by = 'sampleID')

# add in fraction ID column
dat[, fraction := sub('(w0|DIM)\\.(.*)\\.(\\d?)', '\\3', SampleID)
    ][, fraction := as.numeric(fraction)]

# Adjust for spin artifact in the C treatment
# adjustment calculated by taking the difference in median of "good" curves and median of "bad" curves
to_adjust <- paste0('W1_', c('PP_4', 'PJ_5', 'MC_6', 'GL_4', 'GL_5', 'PJ_4',
                             'PP_5', 'MC_5', 'MC_23', 'PJ_23', 'GL_24', 'PP_22',
                             'PJ_22', 'MC_22', 'PP_9'))
density_correct <- 0.0123115812168361760115
dat[sampleID %in% to_adjust, Density.g.ml := Density.g.ml - density_correct]
#
rm(density_correct, to_adjust)

# identify some "off-target" lineages to demonstrate filtering capabilities
ex_bad <- c(dat[Kingdom == 'Archaea'][, uniqueN(sampleID), by = .(iso_trt, taxonID)][order(-V1)][, taxonID][1:10],
            tax[Kingdom == 'Eukaryota'][, taxonID],  # only 3 eukaryotes
            tax[Kingdom == 'Bacteria'
                ][, tax_string := paste(Phylum, Class, Order, Family, Genus, Species, sep = ';')
                  ][grepl('mitochond|chloroplast', tax_string, ignore.case = T)][, taxonID][1:15],
            tax[Kingdom == 'Unassigned'][, taxonID][1:5]
)

# limit to just two of the ecosystems to cut down on data size
dat <- dat[ecosystem %in% c('MC', 'GL')]

# identify which organism are highly frequent and abundant in order to limit to these organisms
dat[, tax_freq := uniqueN(sampleID), by = .(taxonID, iso_trt)]
dat <- dat[iso_trt == 'label', tax_freq := uniqueN(sampleID), by = .(taxonID, ecosystem, treatment)]
dat <- dat[order(-tax_freq)]
# frequent taxa in each treatment of the labeled data
label_freq <- dat[iso_trt == 'label'
                  ][, .(freq = uniqueN(sampleID)), by = .(taxonID, ecosystem, treatment)
                    ][, .SD[1:25], by = .(ecosystem, treatment)
                      ][, unique(taxonID)]
# frequent taxa in the light data
light_freq <- dat[iso_trt == 'light'][!duplicated(taxonID)][, taxonID][1:25]
# combine
ex_good <- union(label_freq, light_freq)

# define example qSIP based on these 120 taxa
example_qsip <- dat[taxonID %in% c(ex_bad, ex_good)]

# final cleaning and column ordering
example_qsip <- setcolorder(example_qsip, neworder = c('taxonID', 'sampleID', 'fraction',
                                                       'timepoint', 'isotope', 'iso_trt',
                                                       'ecosystem', 'treatment', 'rep',
                                                       'Density.g.ml', 'avg_16S_g_soil', 'seq_abund',
                                                       'Kingdom', 'Phylum', 'Class', 'Order',
                                                       'Family', 'Genus'))

example_qsip <- example_qsip[, c('taxonID', 'sampleID', 'fraction', 'timepoint',
                                 'isotope', 'iso_trt', 'ecosystem', 'treatment', 'rep',
                                 'Density.g.ml', 'avg_16S_g_soil', 'seq_abund',
                                 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')]

# change the iso_trt values for the t0 samples to NA
example_qsip[timepoint == 0, iso_trt := NA]

# round up 16S abundances
example_qsip[, avg_16S_g_soil := ceiling(avg_16S_g_soil)]

# set factor levels for iso_trt
example_qsip[, iso_trt := factor(iso_trt, levels = c('light', 'label'))
             ][, treatment := factor(treatment,
                                     levels = c('control', 'C', 'C_N'),
                                     labels = c('Control', 'C', 'CN'))]

# re-order by samples and treatments
example_qsip <- example_qsip[order(timepoint, iso_trt, treatment, ecosystem, rep)]

# relabel ASVs based on Phylum and abundance
asv_labels <- example_qsip[, .(tot_abund = sum(seq_abund)), by = .(taxonID, Kingdom, Phylum, Class, Order, Family, Genus)
                           ][, phy_abund := sum(tot_abund), by = Phylum
                             ][order(-phy_abund)]

phy_order <- asv_labels[order(Kingdom)][!duplicated(Phylum), Phylum]

asv_labels[, Phylum := factor(Phylum, levels = phy_order)]
asv_labels <- asv_labels[order(Phylum, Class, -tot_abund)]

asv_labels[, asv_id := paste0('ASV_', 1:.N)]

# combine with example data
example_qsip <- merge(example_qsip, asv_labels[, .(taxonID, asv_id)], by = 'taxonID', all.x = T, sort = F)

setcolorder(example_qsip, 'asv_id')
example_qsip[, taxonID := NULL]

# output final data into /data folder of package directory to be documented and exported
# NOTE THAT THIS IS A DATA.TABLE OBJECT AND WILL BE SAVED AS AN .RDA FILE
usethis::use_data(example_qsip, overwrite = TRUE)
