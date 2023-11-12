#' Title
#'
#' @param file_dir string, the directory storing the mzML files
#' @param mzmlfiles_pattern string, file format of the raw lcms data
#' @param aligned_tb data frame, aligned feature table
#' @param int_threshold numeric, intensity threshold
#'
#' @return data frame
#' @export
#' @import xcms
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#'
#'
fillGap <- function(file_dir, mzmlfiles_pattern, aligned_tb, int_threshold = 0 ){
  align_tb <- aligned_tb
  mzml <- list.files(pattern = mzmlfiles_pattern, path = file_dir)
  start_i <- ncol(align_tb) - length(mzml)
  missing_eic_dir <- paste0(file_dir,"/missing_eic_int_threshold_", int_threshold)
  dir.create(missing_eic_dir)
  recovered_num <- 0
  for(i in 1:length(mzml)){
    index <- which(is.na(align_tb[,start_i+i])  == TRUE)
    lcmsraw <- xcms::xcmsRaw(paste0( file_dir,"/",mzml[i]))

    # j iterate all missing value
    for (j in index){
      t_mz <- align_tb$mz[j]
      t_rt <- align_tb$rt[j]
      eic <- EIC_matrix(lcmsraw, t_mz, t_rt/60)
      eic[,2] <- peak_smooth(eic[,2])
      if(int_threshold == 0){
        int_Thre <- int_threshold +1
      }else{
        int_Thre <- int_threshold
      }
      if(max(eic[,2]) >= int_Thre ){
        recovered_num <- recovered_num + 1
        png(paste0(missing_eic_dir,"/",i,"th_sample_",
                   round(align_tb$mz[j],4),"_",round(align_tb$rt[j]/60),"min.png"))
        plot(eic[,1]/60, eic[,2], type="l",
             xlab="rt, min", ylab="int",
             main = paste0(i,"th_sample_",round(align_tb$mz[j],4),"_",
                           round(align_tb$rt[j]/60),"min_int",round(max(eic[,2]))))
        dev.off()
        align_tb[j,start_i+i] <- max(round(eic[,2]) + 0.01)
      }
    }
  }
  write.csv(align_tb,
            paste0(file_dir,"/",recovered_num,"_gap_filling_with_int_threshold_",int_threshold ,".csv"),
            row.names = FALSE)
  return(align_tb)
}
