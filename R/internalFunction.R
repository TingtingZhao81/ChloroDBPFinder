
#' Title: load MS1 information
#'
#' @param mzMLdirectory string, file directory
#' @param mzMLfile string, file full name .mzML
#'
#' @return list
#' @export
#' @import MSnbase
#'
#'
loadRawMS1 <- function(mzMLdirectory, mzMLfile){
  rawMS1 <- MSnbase::readMSData(paste0(mzMLdirectory,"/",mzMLfile), mode = "onDisk",  msLevel. = 1)
  return(rawMS1)
}



#' Title: record the retention time of MS1 scan
#'
#' @param rawMS1 MSnbase object, list of MS1 spectra
#'
#' @return vector, retention time of all ms1 spectra
#' @export
#' @import MSnbase
#' @import doParallel
#' @import foreach
#' @import parallel
#'
rtMS1scan <- function(rawMS1){
  #internal function

  # Calculate the number of cores
  no_cores <- parallel::detectCores() - 2
  # Initiate cluster
  doParallel::registerDoParallel(no_cores)
  k <- 1

  rt_v <- foreach::foreach( k=1:length(rawMS1), .combine = cbind, .packages = "MSnbase" ) %dopar% {
    return(rawMS1[[k]]@rt)
  }

  rt_v_MS1 <- as.numeric(as.vector(rt_v[1,]))
  #clean up the cluster
  doParallel::stopImplicitCluster()
  return(rt_v_MS1)
}



#' Title: create folder to save EIC and MS1 spectra
#'
#' @param mzMLdirectory string, file directory
#' @param mzMLfile string, file full name .mzML
#'
#' @return string, folder name
#'
#'
createMS1folder <- function(mzMLdirectory, mzMLfile){
  parent_dr <- mzMLdirectory
  if(length(grep('.mzML', mzMLfile) !=0)){
    sample_name <- strsplit(mzMLfile, ".mzML")[[1]]
  }
  if(length(grep('.mzXML', mzMLfile)) != 0){
    sample_name <- strsplit(mzMLfile, ".mzXML")[[1]]
  }
  ms1_graph_dr <- paste0(parent_dr,"/MS1_",sample_name)
  dir.create(ms1_graph_dr)
  return(ms1_graph_dr)
}
