#' Title: prepare raw data for EIC
#'
#' @param mzMLdirectory string, directory of the raw LCMS data
#' @param mzMLfile  string, file name
#'
#' @return eicRawlcms objected extracted by xcmsRaw
#' @export
#' @import xcms
#'
#'
eicRawlcms <- function(mzMLdirectory, mzMLfile) {
  return(xcms::xcmsRaw(paste0(mzMLdirectory,"/",mzMLfile)))
}


#' Title: extract EIC
#'
#' @param rawlcms xcms object processed by xcmsRaw
#' @param mz numeric, mass
#' @param rt numeric, retention time in minutes
#' @param mz_tol numeric, mass tolerance to extract EIC
#' @param rt_tol numeric, retention time window for EIC
#'
#' @return matrix, retention time column and intensity column
#' @export
#' @import xcms
#'
EIC_matrix <- function(rawlcms, mz, rt, mz_tol=0.01, rt_tol=30){
  mzRange <- c(mz-mz_tol, mz + mz_tol)
  rtRange <- c(rt*60-rt_tol, rt*60+rt_tol)
  rawEIC <- xcms::rawEIC(rawlcms, mzrange= mzRange, rtrange = rtRange)
  eic_matrix <- cbind(rawlcms@scantime[rawEIC$scan],
                      rawEIC$intensity)
  return(eic_matrix)
}

### function 7: peak smoothing function ###
#' Title: smooth the chromatographic peaks
#'
#' @param x vector, a list of intensity information of the peak
#' @param level numeric, default 2
#'
#' @return vector, a list of intensity
#' @export
#'
peak_smooth <- function(x,level=2){
  n <- level
  if(length(x) < 2*n){
    return(x)
  } else if(length(unique(x))==1){
    return(x)
  } else{
    y <- vector(length=length(x))
    for(i in 1:n){
      y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
    }
    for(i in (n+1):(length(y)-n)){
      y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
    }
    for(i in (length(y)-n+1):length(y)){
      y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
    }
    return(y)
  }
}
