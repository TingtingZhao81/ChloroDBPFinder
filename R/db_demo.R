#' Demon of database to match the tandem ms2 spectra
#'
#' Contains 8 qualities of the chemicals
#'
#' @format A data frame with 1 row and 8 variables
#' \describe{
#'      \item{Name}{chemical names}
#'      \item{Precursor_type}{adduct form of precursors}
#'      \item{Collision_energy}{collision energy for fragmentation }
#'      \item{Precursor_mz}{mz value of the precursors}
#'      \item{Formula}{molecular formula}
#'      \item{Ion_mode}{ion mode}
#'      \item{MS2mz}{ms2 fragment mz}
#'      \item{MS2int}{ms2 fragment intensity}
#'      }
#' @source {Created in-house to serve as an example}
#'
#' @examples
#' data(db_demo) # Lazy loading. Data becomes visible as soon as
"db_demo"

