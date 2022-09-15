#' Rhesus Macaque clonal tracking dataset
#'
#' A dataset containing clonal tracking cell counts from a Rhesus Macaque study.
#'
#' @format A list containing clonal tracking data for each animal (ZH33, ZH17, ZG66).
#' Each clonal tracking dataset is a 3-dimensional array whose dimensions identify
#' \describe{
#'   \item{1}{time, in months}
#'   \item{2}{cell types: T, B, NK, Macrophages(M) and Granulocytes(G)}
#'   \item{3}{unique barcodes (clones)}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3979461/bin/NIHMS567927-supplement-02.xlsx}
"Y_RM"

#' Clonal tracking data from clinical trials
#'
#' A dataset containing clonal tracking cell counts from three different clinical trials.
#'
#' @format A list containing the clonal tracking data for each clinical trial (WAS, \eqn{\beta 0 \beta E}{b0bE}, \eqn{\beta S \beta S}{bSbS}).
#' Each clonal tracking dataset is a 3-dimensional array whose dimensions identify
#' \describe{
#'   \item{1}{time, in months}
#'   \item{2}{cell types: T, B, NK, Macrophages(M) and Granulocytes(G)}
#'   \item{3}{unique barcodes (clones)}
#' }
#' @source \url{https://github.com/BushmanLab/HSC_diversity}
"Y_CT"
