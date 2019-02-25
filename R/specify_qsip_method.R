#' Specify qSIP-related data to phylosip object
#'
#' Creates or updates a phylosip object with experimental data and/or filtering criteria
#'
# @param formula Formula specifying density, abundance, and sample ID data in the following order: \code{density + abund ~ rep_id + rep_group + timepoint}
#' @param data Data as a \code{phyloseq} object
#' @param density Single length character matching to variable in \code{data} describing sample densities in UNITS
#' @param abund Single length character matching to variable in \code{data} describing sample 16S gene abundance
#' @param rep_id Single length character matching to variable in \code{data} describing sample/replicate names to summarize across fractions.
#'   Required for calcualtion of weighted average DNA densities (\emph{per sample})
#' @param rep_num Single length character matching to variable in \code{data} matching replicates from the same sample origin across treatments
#'   and any possible replicate groupings specified by \code{rep_group}.
#'   Required if calculating \emph{per-capita} population rates of change where differences in abundance are only considered between sample from
#'   the same origin
#' @param rep_group Single length character matching to variable in \code{data} used to summarize values across groups of replicates.
#'   Required for calculation of atom excess fraction
#' @param iso Single length character vector describing which isotope was used, choices of \code{18O} or \code{13C} describing either
#'   \eqn{^18^O}-labeled water or \eqn{^13^C}-labeled carbon inputs, respectively.
#'   Required for calculation of molecular weights and atom excess fraction
#' @param iso_trt Single length character matching to variable in \code{data} used to distinguish between samples with and without isotope added.
#'   The data this value references should be a facotr with "light" treatments as the first level.
#' @param timepoint Single length character matching to variable in \code{data} used to distinguish between samples after certain incubation
#'   timepoints. The data that this value references should be a factor with time 0 values as the first level.
#'   Required for calculation of population birth and death rates
#' @param filter Optional data frame of filtering options to apply to taxa based on sample and fraction-level frequency
#' @param ... Arguments to change variables of existing \code{phylosip} object
#'
#' @details Some details
#'
#' @return \code{specify_qsip} returns a \code{phyloseq}-like object of class \code{phylosip}.
#'   The \code{phylosip} object is an extension of \code{phyloseq} in that it contains an additional formal class,
#'   the \code{qsip} class. Objects of \code{qsip} class contain 7 slots and a \code{.Data} slot which is a list of
#'   all qSIP-related calculations (weighted average densities, molecular weights, and/or atom excess fraction) as well
#'   as any corresponding confidence intervals.
#'
#' @examples
#'  # Create a new phylosip object from a phyloseq object
#'
#'  # Update an existing phylosip object
#'
#' @export

setGeneric('specify_qsip',
           valueClass=c('phyloseq', 'phylosip'),
           function(data,
                    density='character', abund='character', rep_id='character',
                    rep_group='character', iso='character', iso_trt='character',
                    timepoint='character', filter_levels='data.frame', filter='character', ...) {
             standardGeneric('specify_qsip')
           })

setMethod('specify_qsip',
          signature(data='phyloseq'),
          function(data, density=character(), abund=character(), rep_id=character(), rep_num=character(), rep_group=character(),
                   iso=c('18O', '13C'), iso_trt=character(), timepoint=character(), filter_levels=NULL, filter=character()) {
            if(missing(data)) stop('Must supply data as phyloseq object', call.=F)
            #            if(is.null(formula)) { # User input good (w/o formula)?
            if(missing(density) || missing(abund)) stop('Must supply both sample densities and abundances', call.=F)
            # are supplied terms of single character length?
            matching_args <- as.list(environment()) # return function arguments supplied by user
            matching_args <- matching_args[!names(matching_args) %in% c('data', 'filter')]
            input_length <- sapply(matching_args, length)
            if(any(input_length > 1)) {
              culprits <- names(matching_args)[which(multiple_input)]
              stop('Arguments ', paste0(culprits, collapse=', '), 'are not of length 1', call.=F)
            }
            # do supplied terms match existing data names?
            matching_args <- matching_args[!names(matching_args) %in% 'iso'] # ignore iso from this
            matching_args <- matching_args[input_length==1] # only consider valid inputs
            matching_args <- sapply(matching_args, match, data@sam_data@names)
            if(anyNA(matching_args)) {
              culprits <- names(matching_args)[which(is.na(matching_args))] # who doesn't match
              stop('Arguments ', paste0(culprits, collapse=', '), ' do not mach those in @sam_data object', call.=F)
            }
            #            } else {  # User input good (w formula)?
            #              formula_terms <- eval_formula(formula)
            #              if(length(formula_terms$response < 2 || formula_terms$response > 2)) {
            #                stop('Must specify both density and abundance measurements as response variables, respectively', call.=F)
            #              }
            #              density_name_match <- match(formula_terms$response[1], data@sam_data@names)
            #              abund_name_match <- match(formula_terms$response[2], data@sam_data@names)
            #            }
            #              if(length(formula_terms$grouping)>=1) rep_id <- formula_terms$response[1]
            #              if(length(formula_terms$grouping)>=2) rep_group <- formula_terms$response[2]
            #              if(length(formula_terms$grouping)>=3) timepoint <- formula_terms$response[3]
            iso <- match.arg(toupper(iso), c('18O', '13C'))
            data <- as(data, 'phylosip')
            if(!missing(density)) data@qsip@density <- density
            if(!missing(abund)) data@qsip@abund <-  abund
            if(!missing(rep_id)) data@qsip@rep_id <- rep_id
            if(!missing(rep_num)) data@qsip@rep_num <- rep_num
            if(!missing(rep_group)) data@qsip@rep_group <- rep_group
            if(!missing(timepoint)) data@qsip@timepoint <- timepoint
            if(!missing(iso)) data@qsip@iso <- iso
            if(!missing(iso_trt)) data@qsip@iso_trt <- iso_trt
            if(!missing(filter)) data@qsip@filter <- filter
            if(!missing(filter_levels)) data@qsip@filter_levels <- filter_levels
            return(data)
          })

# You can update your qsip object when you call calculation functions by
# simply calling specify_qsip(data, ...) in the calculation function code
setMethod('specify_qsip',
          signature(data='phylosip'),
          function(data, ...) {
            if(!missing(density)) data@qsip@density <- density
            if(!missing(abund)) data@qsip@abund <-  abund
            if(!missing(rep_id)) data@qsip@rep_id <- rep_id
            if(!missing(rep_num)) data@qsip@rep_num <- rep_num
            if(!missing(rep_group)) data@qsip@rep_group <- rep_group
            if(!missing(timepoint)) data@qsip@timepoint <- timepoint
            if(!missing(iso)) data@qsip@iso <- iso
            if(!missing(iso_trt)) data@qsip@iso_trt <- iso_trt
            if(!missing(filter)) data@qsip@filter <- filter
            if(!missing(filter_levels)) data@qsip@filter_levels <- filter_levels
            return(data)
          })
