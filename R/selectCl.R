#' Title: Identify the Cl peaks in a feature table
#'
#' @param mzMLdirectory string, the directory containing the raw mzML/mzXML file
#' @param mzMLfile string, name of the raw mzML/mzXML file
#' @param original_ft data frame, feature table
#' @param binary_model pre-trained binary model to predict the existence of the Cl
#' @param multi_model pre-trained multiclass model to predict the number of Cl
#' @param ms1_spetra_rt_tol numeric, second unit, retention time window to extract the ms1 spectra
#' @param ms1_spectra_mass_tol numeric, Da unit, m/z tolerance to extract the ms1 spectra
#' @param iso_mass_diff_4 numeric, M4 - M0 m/z difference
#' @param iso_mass_diff_3 numeric, M3 - M0 m/z difference
#' @param iso_mass_diff_2 numeric, M2 - M0 m/z difference
#' @param iso_mass_diff_1 numeric, M1 - M0 m/z difference
#'
#' @return data frame with Cl information
#' @export
#' @import xcms
#' @import MSnbase
#' @import doParallel
#' @import foreach
#' @import randomForest
#'
#'
selectCl <- function(mzMLdirectory, mzMLfile, original_ft,
                     binary_model, multi_model,
                     ms1_spetra_rt_tol = 20, ms1_spectra_mass_tol = 25,
                     iso_mass_diff_4 = 3.994 , iso_mass_diff_3 = 3.000,
                     iso_mass_diff_2 = 1.99705, iso_mass_diff_1= 1.003355){

  rawMS1 <- loadRawMS1(mzMLdirectory, mzMLfile)
  ms1_rt_v <- rtMS1scan(rawMS1)
  eiclcms <- eicRawlcms(mzMLdirectory, mzMLfile)
  ms1_graph_dr <- createMS1folder(mzMLdirectory, mzMLfile)

  features <- paste0(original_ft$featureID,";", original_ft$mz, ";",original_ft$rt)
  # Calculate the number of cores
  no_cores <- parallel::detectCores() - 2
  # Initiate cluster
  doParallel::registerDoParallel(no_cores)
  k <- 1
  prediction_tb <- foreach::foreach(k =1:length(features), .combine = rbind, .packages = c("xcms","MSnbase","randomForest")) %dopar% {
    predictCl(features[k], rawMS1, ms1_rt_v, rawlcms= eiclcms,
              ms1_graph_dr =  ms1_graph_dr, binary_model, multi_model,
              ms1_spetra_rt_tol = ms1_spetra_rt_tol, ms1_spectra_mass_tol = ms1_spectra_mass_tol,
              iso_mass_diff_4 = iso_mass_diff_4 , iso_mass_diff_3 = iso_mass_diff_3,
              iso_mass_diff_2 = iso_mass_diff_2, iso_mass_diff_1= iso_mass_diff_1)
  }
  nrow(prediction_tb)

  #clean up the cluster
  doParallel::stopImplicitCluster()

  prediction_tb <- as.data.frame(prediction_tb)
  colnames(prediction_tb) <-  c("cl", "M0_mz","M1_mz","M2_mz","M3_mz", "M4_mz",
                                "M0_int", "M1_int","M2_int","M3_int","M4_int","cor")
  cl_tb <- cbind(original_ft, prediction_tb)
  return(cl_tb)
}



#' Title: extract and predict individual peak
#'
#' @param featureInfo string,"feature ID", "m/z" and "RT" separated by ";"
#' @param rawMS1 MS1 spectra object extracted using MSnbase
#' @param ms1_rt_v list, retention time list
#' @param rawlcms raw LCMS extracted by xcmsRaw function
#' @param ms1_graph_dr folder to save MS1 and EIC
#' @param binary_model pre-trained binary model to predict the existence of the Cl
#' @param multi_model pre-trained multiclass model to predict the number of Cl
#' @param ms1_spetra_rt_tol numeric, second unit, retention time window to extract the ms1 spectra
#' @param ms1_spectra_mass_tol  numeric, Da unit, m/z tolerance to extract the ms1 spectra
#' @param iso_mass_diff_4 numeric, M4 - M0 m/z difference
#' @param iso_mass_diff_3 numeric, M3 - M0 m/z difference
#' @param iso_mass_diff_2 numeric, M2 - M0 m/z difference
#' @param iso_mass_diff_1 numeric, M1 - M0 m/z difference
#'
#' @return prediction result and MS1 m/z column and RT column
#' @import xcms
#' @import MSnbase
#' @import doParallel
#' @import foreach
#' @import randomForest
#' @importFrom stats cor
#' @importFrom stats sd
#' @importFrom stats predict
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#'
predictCl <- function(featureInfo, rawMS1, ms1_rt_v, rawlcms, ms1_graph_dr,
                      binary_model, multi_model,
                      ms1_spetra_rt_tol = 20, ms1_spectra_mass_tol = 25,
                      iso_mass_diff_4 = 3.994 , iso_mass_diff_3 = 3.000,
                      iso_mass_diff_2 = 1.99705, iso_mass_diff_1= 1.003355){
  # extract information
  featureInfo <- strsplit(featureInfo, ";")[[1]]
  ID <- featureInfo[1]
  t_mz <- as.numeric(featureInfo[2])
  t_rt <- as.numeric(featureInfo[3])
  # prepare table to record
  binary_features <- c( "mz0", "int2_o_int0", "int1_o_int0", "RI2_RI1","mz_1_0","mz_2_0")
  binary_table <- as.data.frame(matrix(ncol=length(binary_features), nrow=1))
  colnames(binary_table) <- binary_features
  multi_features <- c( "mz0","int2_o_int0", "int1_o_int0", "RI2_RI1","mz_1_0","mz_2_0",
                       "mz_3_0", "mz_4_0","int3_o_int0", "int4_o_int0", "RI4_RI3")
  multi_table <- as.data.frame(matrix(ncol=length(multi_features), nrow=1))
  colnames(multi_table) <- multi_features

  # 3.1 locate MS1 in the raw data #
  eic <- EIC_matrix(rawlcms, t_mz, t_rt/60, rt_tol = ms1_spetra_rt_tol)
  # int_index <- which(eic[,2] >=0)
  eic[,2] <- peak_smooth(eic[,2])
  rt_v <- eic[,1]
  int_v <- eic[,2]
  normalized_range <- ( max(eic[,2]) - min(eic[,2]))/mean(eic[,2])
  normalized_sd <- ( sd(eic[,2]))/mean(eic[,2])

  ############# evaluate peak shape quality #######
  #### bad peak quality, just stop and return 0 cl
  if(normalized_range <= 0.9 & normalized_sd <= 0.25){
    return(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ,0))
  }
  ##### start from the apex, find the start and stop point ####
  max_index <- which.max(int_v)
  # some extra case that the boundary larger than the apex,
  # correct the max index
  index <- which.min(abs(rt_v - t_rt))
  if(abs(rt_v[index] - rt_v[max_index]) >= 1/2 * ms1_spetra_rt_tol){
    max_index <- index
  }

  # find the range of rt
  # right
  right_index <- c()
  int_prior <- int_v[max_index]
  if((max_index+1) <= length(int_v)){
    for(r in (max_index+1):length(int_v)){
      if(int_v[r] < 1000 | int_v[r] >= int_prior){
        break
      }
      right_index <- c(right_index, r)
      int_prior <- int_v[r]
    }
  }

  # points(rt_v[right_index]/60, int_v[right_index], col="red")
  left_index <- c()
  int_prior <- int_v[max_index]
  if((max_index-1) >= 1){
    for(l in (max_index-1):1){
      if(int_v[l] < 1000 | int_v[l] >= int_prior){
        break
      }
      left_index <- c(left_index, l)
      int_prior <- int_v[l]
    }
  }

  #points(rt_v[left_index]/60, int_v[left_index], col="blue")
  rt_v <- rt_v[sort(c(left_index, max_index, right_index))]

  int_v <- int_v[sort(c(left_index, max_index, right_index))]
  #points(rt_v, int_v, col="purple")

  rt_index <- which(ms1_rt_v %in% rt_v)
  M0_mz_v <- c()
  M0_int_v <- c()
  M1_mz_v <- c()
  M1_int_v <- c()
  M2_mz_v <- c()
  M2_int_v <- c()
  M3_mz_v <- c()
  M3_int_v <- c()
  M4_mz_v <- c()
  M4_int_v <- c()
  # extract the int and mz information in the located ms1 spectra
  for(m in 1:length(rt_index)){
    ms1_spectra <- rawMS1[[rt_index[m]]]
    # M+0
    M0mz_index <- which(abs(ms1_spectra@mz - t_mz) <= t_mz * ms1_spectra_mass_tol *0.000001)
    if(length(M0mz_index) == 0){
      M0_mz_v <- c(M0_mz_v, t_mz)
      M0_int_v <- c(M0_int_v , 0)
    }else{
      int0_index <- which.max(ms1_spectra@intensity[M0mz_index])
      if(length(int0_index) >= 1){
        int0_index <- int0_index[1]
      }
      M0_index <- M0mz_index[int0_index]
      M0_mz_v <- c(M0_mz_v, ms1_spectra@mz[M0_index])
      M0_int_v <- c(M0_int_v , ms1_spectra@intensity[M0_index])
    }

    # M+1
    M1mz_index <- which(abs(ms1_spectra@mz - t_mz - iso_mass_diff_1) <= (t_mz + iso_mass_diff_1 )* ms1_spectra_mass_tol *0.000001)
    if(length(M1mz_index) == 0){
      M1_mz_v <- c(M1_mz_v, t_mz + iso_mass_diff_1 )
      M1_int_v <- c(M1_int_v , 0)
    }else{
      int1_index <- which.max(ms1_spectra@intensity[M1mz_index])
      if(length(int1_index) >= 1){
        int1_index <- int1_index[1]
      }
      M1_index <- M1mz_index[int1_index]
      M1_mz_v <- c(M1_mz_v, ms1_spectra@mz[M1_index])
      M1_int_v <- c(M1_int_v , ms1_spectra@intensity[M1_index])
    }

    # M+2
    M2mz_index <- which(abs(ms1_spectra@mz - t_mz - iso_mass_diff_2) <= (t_mz  + iso_mass_diff_2)* ms1_spectra_mass_tol *0.000001)
    if(length(M2mz_index) == 0){
      M2_mz_v <- c(M2_mz_v, t_mz + iso_mass_diff_2 )
      M2_int_v <- c(M2_int_v , 0)
    }else{
      int2_index <- which.max(ms1_spectra@intensity[M2mz_index])
      if(length(int2_index) >= 1){
        int2_index <- int2_index[1]
      }
      M2_index <- M2mz_index[int2_index]
      M2_mz_v <- c(M2_mz_v, ms1_spectra@mz[M2_index])
      M2_int_v <- c(M2_int_v , ms1_spectra@intensity[M2_index])
    }

    # M+3
    M3mz_index <- which(abs(ms1_spectra@mz - t_mz - iso_mass_diff_3) <= (t_mz  + iso_mass_diff_3) * ms1_spectra_mass_tol *0.000001)
    if(length(M3mz_index) == 0){
      M3_mz_v <- c(M3_mz_v, t_mz + iso_mass_diff_3 )
      M3_int_v <- c(M3_int_v , 0)
    }else{
      int3_index <- which.max(ms1_spectra@intensity[M3mz_index])
      if(length(int3_index) >= 1){
        int3_index <- int3_index[1]
      }
      M3_index <- M3mz_index[int3_index]
      M3_mz_v <- c(M3_mz_v, ms1_spectra@mz[M3_index])
      M3_int_v <- c(M3_int_v , ms1_spectra@intensity[M3_index])
    }

    # M+4
    M4mz_index <- which(abs(ms1_spectra@mz - t_mz - iso_mass_diff_4) <= (t_mz  + iso_mass_diff_4) * ms1_spectra_mass_tol *0.000001)
    if(length(M4mz_index) == 0){
      M4_mz_v <- c(M4_mz_v, t_mz + iso_mass_diff_4 )
      M4_int_v <- c(M4_int_v , 0)
    }else{
      int4_index <- which.max(ms1_spectra@intensity[M4mz_index])
      if(length(int4_index) >= 1){
        int4_index <- int4_index[1]
      }
      M4_index <- M4mz_index[int4_index]
      M4_mz_v <- c(M4_mz_v, ms1_spectra@mz[M4_index])
      M4_int_v <- c(M4_int_v , ms1_spectra@intensity[M4_index])
    }

  }
  M0_int_v <- peak_smooth(M0_int_v)
  M1_int_v <- peak_smooth(M1_int_v)
  M2_int_v <- peak_smooth(M2_int_v)
  M3_int_v <- peak_smooth(M3_int_v)
  M4_int_v <- peak_smooth(M4_int_v)

  # int of M+0 should not be lower than 1/5 apex intensity
  M0int_index <- which(M0_int_v >= 1/5 * max(M0_int_v))
  int0 <- mean(M0_int_v[M0int_index])

  int1 <- mean(M1_int_v[M0int_index])
  # correlation is based on the whole range of EIC #
  cor1_o_0 <- round(cor(M0_int_v, M1_int_v),2)
  int1_o_0 <- int1/int0

  int2 <- mean(M2_int_v[M0int_index])
  cor2_o_0 <- round(cor(M0_int_v, M2_int_v),2)
  int2_o_0 <- int2/int0

  int3 <- mean(M3_int_v[M0int_index])
  cor3_o_0 <- round(cor(M0_int_v, M3_int_v),2)
  int3_o_0 <- int3/int0

  int4 <- mean(M4_int_v[M0int_index])
  cor4_o_0 <- round(cor(M0_int_v, M4_int_v),2)
  int4_o_0 <- int4/int0

  ms1_mz <- c(mean(M0_mz_v[M0int_index]), mean(M1_mz_v[M0int_index]),
              mean(M2_mz_v[M0int_index]), mean(M3_mz_v[M0int_index]), mean(M4_mz_v[M0int_index]))
  ms1_int <-  c(int0, int1, int2, int3, int4 )

  cor_1_4 <- paste0("M+1:",cor1_o_0, "_M+2:", cor2_o_0,
                    "_M+3:", cor3_o_0,"_M+4:", cor4_o_0)

  result_list <- c(round(ms1_mz,4), round(ms1_int,2), cor_1_4)
  # for some extra intensity value, skip, note as unqualified
  if(length(which(c(int1_o_0, int2_o_0, int3_o_0, int4_o_0) >= 2)) != 0){
    return(c(0, result_list))
  }
  # correlation of M+2 with M+0 smaller than -0.1, skip, recorded as non-cl
  if(!is.na(cor2_o_0)){
    if(cor2_o_0 <= -0.1){
      return(c(0, result_list))
    }
  }
  # make predictions by the first model
  binary_table[1,] <- as.numeric( c(t_mz, int2_o_0, int1_o_0, int2_o_0 - int1_o_0,
                                    ms1_mz[2] - ms1_mz[1], ms1_mz[3] - ms1_mz[1]))
  pre1 <- predict(binary_model, binary_table )
  # non-cl, skip to next feature

  if( pre1 == 0 ){
    return(c(0, result_list))
  }

  # cl compounds, do the second prediction
  # make predictions using the second model
  multi_table[1,] <- as.numeric( c(t_mz, int2_o_0, int1_o_0, int2_o_0 - int1_o_0,
                                   ms1_mz[2] - ms1_mz[1], ms1_mz[3] - ms1_mz[1],
                                   ms1_mz[4] - ms1_mz[1], ms1_mz[5] - ms1_mz[1],
                                   int3_o_0, int4_o_0, int4_o_0 -int3_o_0 ))

  pre2 <- predict(multi_model,multi_table )

  # only output the EIC and MS1 for Cl-compounds

  png(paste0(ms1_graph_dr,"/feature",ID,"_Cl_",pre2,
             "_mz",round(t_mz,4),
             "_rt",round(t_rt/60,0),".png"))
  ymax <- max(c(eic[,2],M1_int_v ,M2_int_v, M3_int_v,M4_int_v))
  plot(eic[,1]/60, eic[,2],
       ylim=c(0, ymax),
       xlim=c(t_rt/60 -0.3, t_rt/60 +0.3),
       type="l",
       xlab="rt, min",
       ylab="intensity",
       main=paste("feature", ID,"Cl",pre2,"mz",round(t_mz,4), "rt",round(t_rt/60,0),
                  "ratios",round(int1_o_0,3),round(int2_o_0,3),
                  round(int3_o_0,3), round(int4_o_0,3)))
  points(eic[,1]/60, eic[,2])
  # M+0
  lines(rt_v[M0int_index]/60, M0_int_v[M0int_index], type = "l",col="purple")
  points(rt_v/60, M0_int_v, type = "l",col="purple")
  # M+1
  lines(rt_v[M0int_index]/60, M1_int_v[M0int_index], type="l", lty=1 ,col="red")
  points(rt_v/60, M1_int_v, col="red")
  # M+2
  lines(rt_v[M0int_index]/60, M2_int_v[M0int_index], type="l", lty=1 ,col="blue")
  points(rt_v/60, M2_int_v, col="blue")
  # M+3
  lines(rt_v[M0int_index]/60, M3_int_v[M0int_index], type="l", lty=1 ,col="#009999")
  points(rt_v/60, M3_int_v, col="#009999")
  # M+4
  lines(rt_v[M0int_index]/60, M4_int_v[M0int_index], type="l", lty=1 ,col="gray33")
  points(rt_v/60, M4_int_v, col="gray33")
  legend("topright", legend=c("M+0",paste0("M+1, ", round(int1/int0,2)),
                              paste0("M+2, ", round(int2/int0,2)),
                              paste0("M+3, ", round(int3/int0,2)),
                              paste0("M+4, ", round(int4/int0,2))),
         fill= c("purple", "red" ,"blue", "#009999", "gray33"))
  dev.off()



  # output MS1 spectra

  png(paste0(ms1_graph_dr,"/feature",ID,"_Cl_",pre2,"_mz_",round(t_mz,4), "_rt_",round(t_rt/60,0),".png"))
  plot(ms1_mz,
       ms1_int,
       type="h", xlab="m/z", ylab="intensity",
       xlim=c(t_mz-1, t_mz+5),
       ylim=c(0,max(ms1_int)*1.1),
       main = paste("feature",ID,"Cl",pre2,"mz",round(t_mz,4),"rt",round(t_rt/60,3),
                    "ratios",round(int1_o_0,3),round(int2_o_0,3),
                    round(int3_o_0,3), round(int4_o_0,3)))
  text(ms1_mz, ms1_int,labels = as.character(round(ms1_mz,4)),
       pos = 3,cex=0.8)
  dev.off()

  return(c(pre2, result_list))
}
