#' Plot density curves
#'
#' Produces plot of density curves for data frame, phyloseq, or phylosip class objects using ggplot
#'
#' @param data Data as a \code{phylosip} object
#' @param density Single length character matching to variable in \code{data} describing sample densities
#' @param abund Single length character matching to variable in \code{data} describing sample 16S gene abundance
#' @param iso_trt Single length character matching to variable in \code{data} used to distinguish between samples with and without isotope added.
#'   The data this value references should be a factor with "light" treatments as the first level.
#' @param rep_id Single length character matching to variable in \code{data} describing sample/replicate names to summarize across fractions
#' @param rep_group Single length character matching to variable in \code{data} used to summarize values across groups of replicates

#' @param timepoint Single length character matching to variable in \code{data} used to distinguish between samples after certain incubation timepoints

#' @details Only the \code{data} argument is required if a \code{phylosip} object is provided. If \code{data} is of class \code{data.frame} or
#'   \code{phyloseq}, then the \code{density}, \code{abund}, and \code{iso_trt} arguments are required to plot density curves. Other arguments are
#'   required if replicates are to be grouped.
#'   Rows with values of \code{NA} in any column will be omitted.
#'
#' @return \code{plot_curve} returns a ggplot object that can be customized and added to using the common ggplot syntax.
#'
#' @examples
#'  # Load in example data
#'
#'  # plot curve by specifying columns
#'
#'  # specify qsip
#'
#'  # plot curve
#'
#' @export

setGeneric('plot_curve',
           valueClass=c('data.frame', 'phyloseq', 'phylosip'),
           function(data,
                    density='character', abund='character', iso_trt='character', rep_id='character',
                    rep_group='character', timepoint='character', vertical='logical', ...) {
             standardGeneric('plot_curve')
           })

setMethod('plot_curve',
          signature(data='phylosip'),
          function(data, vertical=FALSE) {
            if(missing(data)) stop('Must supply data', call.=F)
            if(length(data@qsip@iso_trt)==0) stop('Must supply isotope specifications', call.=F)
            # create data frame including all specified and relevant columns
            return_columns <- c(data@qsip@density,
                                data@qsip@abund,
                                data@qsip@iso_trt,
                                ifelse(length(data@qsip@rep_id)!=0, data@qsip@rep_id, ''),
                                ifelse(length(data@qsip@rep_group)!=0, data@qsip@rep_group, ''),
                                ifelse(length(data@qsip@timepoint)!=0, data@qsip@timepoint, ''))
            return_columns <- return_columns[nchar(return_columns) > 0]
            plot_data <- as(data@sam_data[,return_columns], 'data.frame')
            # omit NA rows
            plot_data <- plot_data[complete.cases(plot_data),]
            # create basic plot of abundance by density
            basic_plot <- ggplot(plot_data, aes_string(x=data@qsip@density, y=data@qsip@abund, color=data@qsip@iso_trt)) +
              ylab('Abundance') +
              xlab('Density')
            if(vertical) basic_plot + coord_flip() + scale_x_reverse()
            return(basic_plot)
          })

setMethod('plot_curve',
          signature(data='phyloseq'),
          function(data, density=character(), abund=character(), iso_trt=character(), rep_id=character(),
                   rep_group=character(), timepoint=character(), vertical=FALSE) {
            if(missing(data)) stop('Must supply data', call.=F)
            if(missing(density) || missing(abund) || missing(iso_trt)) stop('Must supply sample densities, abundances, and isotope specifications', call.=F)
            # are supplied terms of single character length?
            matching_args <- as.list(environment()) # return function arguments supplied by user
            matching_args <- matching_args[!names(matching_args) %in% c('data', 'filter')]
            input_length <- sapply(matching_args, length)
            multiple_input <- input_length > 1
            if(any(multiple_input)) {
              culprits <- names(matching_args)[which(multiple_input)]
              stop('Arguments ', paste0(culprits, collapse=', '), 'are not of length 1', call.=F)
            }
            # do supplied terms match existing data names?
            matching_args <- matching_args[!names(matching_args) %in% 'iso'] # ignore iso from this
            matching_args <- matching_args[input_length==1] # only consider valid inputs
            matching_args <- sapply(matching_args, match, data@sam_data@names)
            if(anyNA(matching_args)) {
              culprits <- names(matching_args)[which(is.na(matching_args))] # who doesn't match
              stop('Arguments ', paste0(culprits, collapse=', '), ' do not match column names in object @sam_data slot', call.=F)
            }
            # create data frame including all specified and relevant columns
            return_columns <- c(density, abund, iso_trt, rep_id, rep_group, iso_trt, timepoint)
            return_columns <- return_columns[nchar(return_columns) > 0]
            plot_data <- as(data@sam_data[,return_columns], 'data.frame')
            # omit NA rows
            plot_data <- plot_data[complete.cases(plot_data),]
            # create basic plot of abundance by density
            basic_plot <- ggplot(plot_data, aes_string(x=density, y=abund, color=iso_trt)) +
              ylab('Abundance') +
              xlab('Density')
            if(vertical) basic_plot + coord_flip() + scale_x_reverse()
            return(basic_plot)
          })

setMethod('plot_curve',
          signature(data='data.frame'),
          function(data, density=character(), abund=character(), iso_trt=character(), rep_id=character(),
                   rep_group=character(), timepoint=character(), vertical=FALSE) {
            if(missing(data)) stop('Must supply data', call.=F)
            if(missing(density) || missing(abund) || missing(iso_trt)) stop('Must supply sample densities, abundances, and isotope specifications', call.=F)
            # are supplied terms of single character length?
            matching_args <- as.list(environment()) # return function arguments supplied by user
            matching_args <- matching_args[!names(matching_args) %in% c('data', 'filter')]
            input_length <- sapply(matching_args, length)
            multiple_input <- input_length > 1
            if(any(multiple_input)) {
              culprits <- names(matching_args)[which(multiple_input)]
              stop('Arguments ', paste0(culprits, collapse=', '), 'are not of length 1', call.=F)
            }
            # do supplied terms match existing data names?
            matching_args <- matching_args[!names(matching_args) %in% 'iso'] # ignore iso from this
            matching_args <- matching_args[input_length==1] # only consider valid inputs
            matching_args <- sapply(matching_args, match, data@sam_data@names)
            if(anyNA(matching_args)) {
              culprits <- names(matching_args)[which(is.na(matching_args))] # who doesn't match
              stop('Arguments ', paste0(culprits, collapse=', '), ' do not match column names in data', call.=F)
            }
            # create data frame including all specified and relevant columns
            return_columns <- c(density, abund, iso_trt, rep_id, rep_group, iso_trt, timepoint)
            return_columns <- return_columns[nchar(return_columns) > 0]
            plot_data <- data[,return_columns]
            # omit NA rows
            plot_data <- plot_data[complete.cases(plot_data),]
            # create basic plot of abundance by density
            basic_plot <- ggplot(plot_data, aes_string(x=density, y=abund, color=iso_trt)) +
              ylab('Abundance') +
              xlab('Density')
            if(vertical) basic_plot + coord_flip() + scale_x_reverse()
            return(basic_plot)
          })
