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

# limit to just two of the ecosystems to cut down on data size

# identify which organism are highly frequent and abundant in order to limit to these organisms


# output final data into /data folder of package directory to be documented and exported
usethis::use_data(example_qsip, overwrite = TRUE)
