qsip Package
================
Bram Stone

-   [Background](#background)
    -   [Terminology](#terminology)
-   [Installation](#installation)
-   [Data preparation](#data-preparation)
    -   [Creation of phyloseq object](#creation-of-phyloseq-object)
-   [Creation of phylosip object](#creation-of-phylosip-object)
    -   [Plot enrichment curves](#plot-enrichment-curves)
-   [Specify frequency filtering criteria](#specify-frequency-filtering-criteria)
-   [Calculate enrichment and growth](#calculate-enrichment-and-growth)
-   [Extract data frames from phylosip object](#extract-data-frames-from-phylosip-object)

Background
----------

Quantitative stable isotope probing (qSIP) is the combination of stable isotope probing – a foundational technique in the study of ecosystems – with microbial targeted amplicon sequencing data.

### Terminology

Conducting a qSIP experiment is, in many ways, an exercise in data organization. While a simple amplicon sequencing study usually generates a single sequencing sample for every individual sample being surveyed, qSIP experiments break each field or experimental sample into usually more than a dozen fractions which must all be sequenced separately. Because of this, as well as to make this package as easy to use as possible, consistent terminology should be applied to any qSIP experiment.

-   **Replicate** - A single physical sample that has been divided into multiple fractions (usually 10–20). The term *replicate* is used rather than *sample* to avoid confusion during the sequencing preparation process, in which each qSIP fraction must be treated as a separate sample. Furthermore, the term *replicate* makes it clear that these represent the unit of statistical replication and power.
-   **Fraction** - A subsample produced by fractionating a single replicate based on vertical stratification following ultra-centrifugation. *Each fraction must have its own qPCR amplicon abundance measurement as well as a density measurement.*
-   **Feature Table** - A table of abundances or relative frequences for each microbial taxon in each sample. For qSIP analyses, feature tables must include the abundances/frequencies in *each* fraction. This is often called an OTU table, ASV table, or species abundance table.
-   **Labeled** - A replicate that has been treated with a stable isotope (<sup>18</sup>O, <sup>13</sup>C, or <sup>15</sup>N). Often used to differentiate labeled from unlabeled values in output from the qsip package (as well as represented in equations throughout the published literature).
-   **Light** - A replicate that has **not** been treated with a stable isotope. An unlabeled replicate.

Installation
------------

**qsip** is not currently on CRAN. The only way to install **qsip** is through Github. This package depends on [phyloseq](http://joey711.github.io/phyloseq/index.html). See here for issues on [installing phylseq](http://joey711.github.io/phyloseq/install.html).

``` r
# install devtools and BiocManager
install.packages('devtools')
install.packages('BiocManager')

# install phyloseq and qsip
BiocManager::install('phyloseq')
devtools::install_github('bramstone/qsip', ref='separate-sample-values')

library(phyloseq)
library(qsip)
```

Data preparation
----------------

Here, go into how to organize and prepare experimental data, the feature table, and taxonomic data for phyloseq combination. Spend the most time on experimental data.

### Creation of phyloseq object

Microbial abundance tables must be combined with taxonomic and experimental data into a single [phyloseq object](http://joey711.github.io/phyloseq/import-data.html).

Creation of phylosip object
---------------------------

``` r
dat <- specify_qsip(dat,
                    density='Density.g.ml',
                    abund='avg_16S_g_soil',
                    rep_id='sampleID',
                    rep_group='ecosystem',
                    iso='18O',
                    iso_trt='isotope',
                    timepoint='timepoint')

dat
```

### Plot enrichment curves

One of the first diagnostic plots that should be generated from an enrichment experiment is the enrichment curve. Often, this is done even before sequencing efforts begin. The basic idea of the enrichment curve is to plot the change in density of DNA due to the incorporation of stable isotopes.

*With qsip, enrichment curves may be plotted from phyloseq, phylosip, or basic data.frame objects.* Because phylosip objects already have replicates, fractions, densities, and qPCR abundances specified, no other information is necessary. Using phyloseq or data.frame objects requires the specification of this information.

``` r
# plot from phylosip object
plot_curve(dat)

# plot from phyloseq object
plot_curve(phlyo_dat, density='density_column', abund='qPCR_column', iso_trt='isotope_treatment_column')

# plot from data.frame object
plot_curve(exper_dat, density='density_column', abund='qPCR_column', iso_trt='isotope_treatment_column')
```

Specify frequency filtering criteria
------------------------------------

``` r
dat@qsip@filter_levels <- create_filters(2, 5, soft=1)
```

Calculate enrichment and growth
-------------------------------

``` r
# atom percent enrichment (APE)
dat <- calc_excess(dat,
                   separate_label=T,
                   filter=T,
                   correction=T)

# per-capita growth
dat <- calc_pop(dat,
                separate_label=T,
                filter=T,
                correction=T,
                growth_model='exponential')
```

Extract data frames from phylosip object
----------------------------------------

``` r
# extract data related to atom percent "excess" and measures from "label"-ed samples
qsip_enrich <- qsmelt(dat, include='excess|label', regex=T)
head(qsip_enrich)

# extract data related to birth and death
qsip_growth <- qsmelt(dat, include='birth|death', abundance=T, regex=T)
head(qsip_growth)
```
