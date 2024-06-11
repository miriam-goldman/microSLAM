#' The microSLAM package
#'
#' microSLAM is a library and R package for microbiome GLMM models
#'
#' @rdname microSLAM
#' @name microSLAM-package
#' @keywords glmm
#' @aliases microSLAM-package microSLAM
#' @docType _PACKAGE
NULL

## usethis namespace: start
#' @importFrom data.table fread setindexv
#' @importFrom RcppXPtrUtils cppXPtr checkXPtr
#' @import RcppArmadillo
#' @importFrom logr put
#' @importFrom magrittr %>%  %<>% 
#' @importFrom dplyr filter left_join right_join select group_by mutate count rename ungroup summarize inner_join n reframe add_count
#' @import ggplot2
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggExtra ggMarginal
#' @importFrom stringr str_split
#' @importFrom parallelDist parDist
#' @importFrom pheatmap pheatmap
#' @import Matrix
#' @importFrom pROC auc
## usethis namespace: end
NULL
