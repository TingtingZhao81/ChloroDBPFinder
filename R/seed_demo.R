#' Demo of seed list to create the network
#'
#' Contains 9 qualities of the chemicals
#'
#' @format A data frame with 11 row and 9 variables
#' \describe{
#'      \item{featureID}{featureID: precursor ? S1_5 }
#'      \item{label}{seed}
#'      \item{annotation}{chemical name, can be empty}
#'      \item{MF}{molecular formula}
#'      \item{mz}{mz}
#'      \item{rt}{retention time, can be empty}
#'      \item{cl}{number of cl, 0/1/2}
#'      \item{ms2_mz}{ms2 fragment mz separated by ";" }
#'      \item{ms2_int}{ms2 fragment intensity separated by ";"}
#'      }
#' @source {Created in-house to serve as an example}
#'
#' @examples
#' data(seed_demo) # Lazy loading. Data becomes visible as soon as
"seed_demo"

