#' The S4 class for storing qSIP related data in list format
#'
#'  Any resulting data will be stored in the \code{.Data} slot. Dscription of each argument found
#'  in \code{specify_qsip} function documentation.
#'
#' @exportClass qsip
setClass('qsip',
         slots=c(density='character', abund='character',
                 rep_id='character', rep_num='character',
                 rep_group='character',
                 iso='character', iso_trt='character',
                 timepoint='character', filter_levels='data.frameOrNULL',
                 filter='character'),
         contains='list')


#' The S4 class for extending phyloseq objects with qSIP related data
#'
#' @importClassesFrom phyloseq phyloseq
#' @exportClass phylosip
setClass('phylosip',
         contains='phyloseq',
         slots=c(qsip='qsip'))
