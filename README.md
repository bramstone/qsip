qsip Package
================
Bram Stone

- <a href="#legacy-package-code" id="toc-legacy-package-code">Legacy
  package code</a>
- <a href="#background" id="toc-background">Background</a>
  - <a href="#terminology" id="toc-terminology">Terminology</a>
- <a href="#installation" id="toc-installation">Installation</a>
  - <a href="#creation-of-phyloseq-object"
    id="toc-creation-of-phyloseq-object">Creation of phyloseq object</a>
- <a href="#creation-of-phylosip-object"
  id="toc-creation-of-phylosip-object">Creation of phylosip object</a>
  - <a href="#plot-enrichment-curves" id="toc-plot-enrichment-curves">Plot
    enrichment curves</a>
- <a href="#specify-frequency-filtering-criteria"
  id="toc-specify-frequency-filtering-criteria">Specify frequency
  filtering criteria</a>
- <a href="#calculate-enrichment-and-growth"
  id="toc-calculate-enrichment-and-growth">Calculate enrichment and
  growth</a>
- <a href="#extract-data-frames-from-phylosip-object"
  id="toc-extract-data-frames-from-phylosip-object">Extract data frames
  from phylosip object</a>
- <a href="#references" id="toc-references">References</a>

## Legacy package code

This branch is no longer the working (i.e., “main”) branch of the `qsip`
package moving forward. The functions contained here have been saved as
“legacy” code meant to allow access to functionality that has been
previously published on. The motivations for this switch were mainly
that this code could not easily handle different types of experimental
designs. As more groups are conducting qSIP experiments under a variety
of settings and designs, it became more efficient to re-implement the
core calculations using functionality from the `data.table` package as
this package moves forward.

The main branch of the `qsip` package can be found at:
<https://github.com/bramstone/qsip/tree/main>

The following publications have utilized functions from this `legacy`
branch:

- Finley *et al.* 2021. Soil minerals affect taxon-specific bacterial
  growth. *The ISME Journal* 16, 1318–1326. DOI:
  [10.1038/s41396-021-01162-y](https://doi.org/10.1038/s41396-021-01162-y)
- Hayer *et al.*. 2021. Microbes on decomposing litter in streams:
  entering on the leaf or colonizing in the water? *The ISME Journal*
  16, 717–725. DOI:
  [10.1038/s41396-021-01114-6](https://doi.org/10.1038/s41396-021-01114-6)
- Stone *et al.* 2021. Nutrients cause consolidation of soil carbon flux
  to small proportion of bacterial community. *Nature Communications*
  12, 1–9. DOI:
  [10.1038/s41467-021-23676-x](https://doi.org/10.1038/s41467-021-23676-x)
- Metz *et al.* 2023. Microbial growth under drought is confined to
  distinct taxa and modified by potential future climate conditions.
  *Nature Communications* 14. DOI:
  [10.1038/s41467-023-41524-y](https://doi.org/10.1038/s41467-023-41524-y)
- Ruan *et al.* 2023. Elevated temperature and CO<sub>2</sub> strongly
  affect the growth strategies of soil bacteria. *Nature
  Communications* 14. DOI:
  [10.1038/s41467-023-36086-y](https://doi.org/10.1038/s41467-023-36086-y)

## Background

**Quantitative stable isotope probing (qSIP)** is the combination of
stable isotope probing – a foundational technique in the study of
ecosystems – with targeted amplicon sequencing data of microbial
communities. In conventional stable isotope probing (SIP) experiments,
identification of the amount of isotopic enrichment of nucleic acids was
done qualitatively – through visual identification and categorization of
nucleic acids into either “heavy” or “light” regions. The qSIP approach
is to divide a single sample into many different fractions (without *a
priori* categorization) along a gradient of increasing densities, and to
estimate the shift in density of an individual microbial taxon’s nucleic
acids based on it’s abundance across the many fractions (Hungate *et
al.* 2015). Further details on the qSIP methodology may be found in
Purcell *et al.* (2019). The core calculation produces estimates of
every microbial taxon’s proportion of enrichment which are often used as
stand-ins for growth. However,`qsip` can also estimate population
*per-capita* rates of growth and death (Koch *et al.* 2018).

The `"legacy"` `qsip` codebase is built on the data structures of the
[`phyloseq` package](http://joey711.github.io/phyloseq/index.html).
`phyloseq` uses S4 object classification to organize different aspects
of microbial sequencing data into a single object (see Data Preparation,
below). The purpose of `phyloseq` is to minimize the code necessary to
perform common operations on microbial community datasets such as
filtering out, or retaining, certain taxonomic lineages or groups of
samples.

The `qsip` legacy code supports stable isotope experiments using
<sup>18</sup>O, <sup>13</sup>C, and <sup>15</sup>N (Morrissey *et al.*
2018). It is also agnostic towards the amplicon being sequenced. It will
support 16S, 18S, ITS, or other amplicon sequencing data.

### Terminology

Conducting a qSIP experiment is, in many ways, an exercise in data
organization. In an amplicon sequencing study, a single sequencing
sample is usually produced from DNA extracted from a single point of
collection (unless samples are pooled, in which case many points of
collection yield one sequencing sample). In a qSIP experiment, DNA from
each sample is divided into usually more than a dozen fractions which
must all be sequenced separately. Because of this, as well as to make
this package as easy to use as possible, consistent terminology should
be applied to any qSIP experiment.

- **Replicate** - A single physical sample that has been divided into
  multiple fractions (usually 10–20). The term *replicate* is used
  rather than *sample* to avoid confusion during the sequencing
  preparation process, in which each qSIP fraction must be treated as a
  separate sample. Furthermore, the term *replicate* makes it clear that
  these represent the unit of statistical replication and power.
- **Fraction** - A subsample produced by fractionating a single
  replicate based on vertical stratification following
  ultra-centrifugation. *Each fraction must have its own qPCR amplicon
  abundance measurement as well as a density measurement.*
- **Feature Table** - Often called an OTU table, ASV table, or species
  abundance table. A table of abundances or relative frequences for each
  microbial taxon in each sample. For qSIP analyses, feature tables must
  include the abundances/frequencies in *each* fraction. For qSIP
  analyses, this should not be a binary presence-absence table.
- **Labeled** - A replicate that has been treated with a stable isotope
  (<sup>18</sup>O, <sup>13</sup>C, or <sup>15</sup>N). Often used to
  differentiate labeled from unlabeled values in output from the qsip
  package (as well as represented in equations throughout the published
  literature).
- **Light** - A replicate that has **not** been treated with a stable
  isotope. An unlabeled replicate.
- **Atom excess fraction (AEF)** - The proportion \[0–1\] of an
  organism’s nucleic acids that are enriched by a stable isotope.
  Identical to **atom percent excess (APE)** or **atom excess percent
  (AEP)**.

## Installation

`qsip` is not currently on CRAN. The only way to install `qsip` is
through Github. To utilize functionality from the legacy codebase, make
sure you specify `ref = "legacy"` when calling `install_github`. See
here if you encounter issues [installing
phyloseq](http://joey711.github.io/phyloseq/install.html).

``` r
# install devtools and BiocManager
install.packages('devtools')
install.packages('BiocManager')

# install phyloseq and qsip using the utilities on devtools and BiocManager
BiocManager::install('phyloseq')
devtools::install_github('bramstone/qsip', ref = 'legacy')

library(phyloseq)
library(qsip)
```

### Creation of phyloseq object

Feature tables must be combined with taxonomic and experimental data
into a single [phyloseq
object](http://joey711.github.io/phyloseq/import-data.html).

## Creation of phylosip object

``` r
dat <- specify_qsip(dat,
                    density='Density.g.ml',
                    abund='avg_16S_g_soil',
                    rep_id='sampleID',
                    rep_group='ecosystem',
                    iso='18O',
                    iso_trt='isotope',
                    timepoint='timepoint')

# high-level qSIP-related data will display upon printing
dat
```

### Plot enrichment curves

One of the first diagnostic plots that should be generated from an
enrichment experiment is the enrichment curve. Often, this is done even
before sequencing efforts begin. The basic idea of the enrichment curve
is to plot the change in density of DNA due to the incorporation of
stable isotopes.

**Note:** previous versions of this README demonstrated the `plot_curve`
function. This function was merely a wrapper to put data into a format
for plotting with `ggplot2`. To utilize the code below, extract the
`@sam_data` slot from your `phyloseq` object, then replace the stand-in
column names with the actual columns that correspond to the fraction
densities and abundances (e.g., qPCR amplicon abundances or DNA
concentration).

``` r
library(ggplot2)

# extract metadata
met <- dat@sam_data

# plot density curves
ggplot(met, aes(frac_density_column, frac_abund_column)) +
  geom_line(aes(group = sample_id_column)) +
  geom_point()
```

## Specify frequency filtering criteria

``` r
dat@qsip@filter_levels <- create_filters(2, 5, soft=1)
```

## Calculate enrichment and growth

``` r
# atom excess fraction (AEF) or atom percent excess (APE)
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

# show new data has been added
dat
names(dat@qsip)
```

## Extract data frames from phylosip object

Previous versions of this README instructed users to use the `qsmelt`
function which was intended to be an analog to phyloseq’s `psmelt`
function and where users could easily match their output to important
experimental variables. However, I did not generate the different
methods required to extract qSIP values that were abstracted to
different levels (i.e., WAD values generated per-sample vs. bootstrapped
enrichment or growth rates generated across many samples – both of which
could exist in the same `phylosip` object).

Instead I recommend users to extract the qSIP values themselves and
merge them into their data as they see fit.

``` r
# extract data related to isotopic "excess" and measures from "label"-ed samples
qsip_enrich <- dat@qsip[['atom_excess']]
head(qsip_enrich)

# extract data related to birth and death
qsip_growth <- dat@qsip[['birth_rate']]
qsip_death <- dat@qsip[['death_rate']]
head(qsip_growth)
head(qsip_death)
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Hungate2015" class="csl-entry">

Hungate BA, Mau RL, Schwartz E *et al.* <span
class="nocase">Quantitative Microbial Ecology through Stable Isotope
Probing</span>. *Applied and Environmental Microbiology* 2015;**81**,
DOI: [10.1128/AEM.02280-15](https://doi.org/10.1128/AEM.02280-15).

</div>

<div id="ref-Koch2018" class="csl-entry">

Koch BJ, McHugh TA, Hayer M *et al.* <span class="nocase">Estimating
taxon-specific population dynamics in diverse microbial
communities</span>. *Ecosphere* 2018;**9**, DOI:
[10.1002/ecs2.2090](https://doi.org/10.1002/ecs2.2090).

</div>

<div id="ref-Morrissey2018" class="csl-entry">

Morrissey EM, Mau RL, Schwartz E *et al.* <span class="nocase">Taxonomic
patterns in the nitrogen assimilation of soil prokaryotes</span>.
*Environmental Microbiology* 2018;**20**, DOI:
[10.1111/1462-2920.14051](https://doi.org/10.1111/1462-2920.14051).

</div>

<div id="ref-Purcell2019" class="csl-entry">

Purcell AM, Dijkstra P, Finley B *et al.* <span
class="nocase">Quantitative Stable Isotope Probing with H218O to Measure
Taxon-Specific Microbial Growth</span>. *Methods of Soil Analysis*
2019;**4**, DOI:
[10.2136/msa2018.0083](https://doi.org/10.2136/msa2018.0083).

</div>

</div>
