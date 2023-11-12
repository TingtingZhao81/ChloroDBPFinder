#' Title: extract chromatographic peaks with MS2 spectra
#'
#' @param mzMLdirectory string, directory of the raw LCMS data
#' @param mzMLfile string, file name
#' @param SN numeric, threshold of peak signal-to-noise. default: 20
#' @param ppm numeric, xcms maximum difference in m/z for peaks with overlapping retention time. default: 20
#' @param peakWidth vector, default: 10 ~ 200
#' @param mzdiff numeric, default: 0.005
#' @param noise numeric,  noise level, peaks with intensity smaller than noise are omitted. default: 2000
#' @param rt_min numeric, lower limit of retention time of interest. default: 300 seconds
#' @param rt_max numeric, upper limit of retention time of interest. default: 3000 seconds
#' @param rt_tol numeric, retention time window in the assignment of tandem MS. default: 10 seconds
#'
#' @return peaks with MS2 spectra
#' @import xcms
#' @import MSnbase
#' @export
#'

extractPeak <- function(mzMLdirectory, mzMLfile, SN = 20, ppm= 20,
                        peakWidth = c(10,120), mzdiff= 0.005,
                        noise=2000, rt_min = 300 , rt_max = 3000, rt_tol = 10){
  lc_ms <- xcms::xcmsSet(paste0(mzMLdirectory,"/",mzMLfile),
                         method = "centWave", # centWave method to extract the feature
                         ppm=ppm, # region of interest: 50 ppm HX recommend
                         peakwidth= peakWidth,
                         mzdiff= mzdiff, # minimum difference in m/z for peaks with overlapping retention time
                         snthresh= SN,
                         prefilter= c(3, 500), # prefilter = c(k, I): mass trace are only retained if they contains
                         # at least k peaks with intensity > I
                         integrate=1,
                         noise=noise # optional argument: peaks with intensity smaller than noise are omitted from ROI
  )

  # peak information: m/z, rt, peak area, peak height
  peaks <- as.data.frame(cbind(1,lc_ms@peaks))
  colnames(peaks)[1] <- "featureID"
  peaks <- peaks[order(peaks$maxo, decreasing = T),]
  peaks$featureID <- 1:nrow(peaks)
  SN_i <- which(peaks$sn >= SN)
  peaks <- peaks[SN_i,]
  peaks <- peaks[peaks$rt <= rt_max,]
  peaks <- peaks[peaks$rt >= rt_min,]

  colnames(peaks)[which(colnames(peaks) == "maxo")] <- "Intensity"
  rownames(peaks) <- peaks$featureID

  ### 2. assign MS2 spectra to each feature
  peaks$ms2_mz <- ""
  peaks$ms2_int <- ""
  rawms2 <- MSnbase::readMSData(paste0(mzMLdirectory,"/",mzMLfile), msLevel. = 2)
  pre_v <- c()
  rt2_v <- c()
  for(num in 1:length(rawms2)){
    pre_v <- c(pre_v, rawms2[[num]]@precursorMz)
    rt2_v <- c(rt2_v,rawms2[[num]]@rt)
  }

  for(i in 1:nrow(peaks)){
    t_mz <- peaks$mz[i]
    t_rt <- peaks$rt[i]

    rt_index <- which(abs(rt2_v - t_rt)<= rt_tol)
    mz_index <- which(abs(pre_v[rt_index] - t_mz) <= 0.01)

    if(length(mz_index) == 0){next}
    # if multiple, assign closet RT
    index <- which.min(abs(rt2_v[rt_index[mz_index]] - t_rt))
    final_index <- rt_index[mz_index][index]

    ms2_mz <- rawms2[[final_index]]@mz
    ms2_int <- rawms2[[final_index]]@intensity
    # filter out noise ( <1% )
    int_filter <- which(ms2_int/max(ms2_int) >= 0.01)
    peaks$ms2_mz[i] <- paste0(round( ms2_mz[int_filter], 4), collapse = ";")
    peaks$ms2_int[i] <- paste0(round( ms2_int[int_filter], 4), collapse = ";")
    #plot(ms2_mz[int_filter], ms2_int[int_filter], type="h")
  }

  return(peaks)
}
