### dataClean ###
#' Title
#'
#' @param peaks data frame, all peaks
#' @param chlorine_tb data frame, all Cl-containing peaks
#' @param rawlcms lcms object extracted by xcmsRaw
#' @param rawfile_dir string, file folder storing the raw mzML/mzXML files
#' @param lcmsfile string, mzML/mzXML file
#' @param adducts logic, TRUE or FALSE
#' @param isotopes logic, TRUE or FALSE
#' @param inSourceFrag logic, TRUE or FALSE
#' @param samNum numeric, sample number
#' @param rt_tol numeric, retention time in seconds
#' @param mz_tol numeric, retention time in seconds
#' @param ISFtable data frame, the primary in source fragment table
#'
#' @return data frame with all cl features
#'
#' @export
#'
cleanFeature <- function(peaks, chlorine_tb, rawlcms, rawfile_dir, lcmsfile,
                         adducts = TRUE, isotopes = TRUE, inSourceFrag = TRUE,
                         ISFtable,
                         samNum = 1, rt_tol = 10, mz_tol = 0.01  ){

  if(adducts){
    cl_tb <- findAdduct(chlorine_tb, rawlcms, samNum,rt_tol ,mz_tol)
  }else{
    cl_tb <- cbind(chlorine_tb,0)
    colnames(cl_tb)[ncol(cl_tb)] <- "Adduct"
  }

  if(isotopes){
    cl_tb <- findIsotope(cl_tb, rawlcms, samNum,rt_tol ,mz_tol)
  }else{
    cl_tb <- cbind(cl_tb,0)
    colnames(cl_tb)[ncol(cl_tb)] <- "isotope"
  }

  if(inSourceFrag){
    cl_tb <- findInsourcefrag(peaks, cl_tb, rawfile_dir,lcmsfile,
                              sampleNumber = samNum, ISFtable = ISFtable)
  }else{
    cl_tb <- cbind(cl_tb,0)
    colnames(cl_tb)[ncol(cl_tb)] <- "ISF_level"
  }
  cl_tb$featureID <- paste0("S", samNum, "_",cl_tb$featureID )
  return(cl_tb)
}

#' Title: identify salt adduct ##
#'
#' @param chlorine_tb data frame, Cl peaks summarized in a table
#' @param rawlcms xcms object processed by xcmsRaw
#' @param samNum numeric
#' @param rt_tol numeric, retention time in seconds
#' @param mz_tol numeric, m/z  in Da
#'
#' @return data frame with adduct labelled
#' @import xcms
#' @importFrom stats cor
#'
findAdduct <- function(chlorine_tb, rawlcms, samNum = 1, rt_tol = 10, mz_tol = 0.01 ){
  peaks <- chlorine_tb
  peaks$Adduct <- 0
  # H+: 1.007276
  # Na+:  22.989221
  NaAdduct <- 22.989221- 1.007276
  # K+: 38.963158
  KAdduct <- 38.963158 - 1.007276
  # NH4+: 18.033826
  NH4Adduct <- 18.033826 - 1.007276
  # -H2O:
  H2OAdduct <- 18.010565
  # + CH3OH
  CH3OHAdduct <- 32.026215

  # iterate for each feature
  for(k in 1:nrow(peaks)){
    t_mz <- peaks$mz[k]
    t_rt <- peaks$rt[k]
    t_num <- peaks$cl[k]
    t_eic_matrix <- EIC_matrix(rawlcms, t_mz, t_rt/60)
    t_int_v <- t_eic_matrix[,2]
    t_int_v <- peak_smooth(t_int_v)
    rt_index <- which(abs(peaks$rt - t_rt) <= rt_tol)

    # we believe the K Na NH4 CH3OH adducts should have lower int than the precursor
    int_col <- which(colnames(peaks) == "maxo")
    if(length(int_col) == 0) {
      int_col <- which(colnames(peaks) == "Intensity")
    }
    int_index <- which(peaks[rt_index,int_col] <= peaks[k,int_col])
    rt_index <- rt_index[int_index]

    if(length(rt_index) == 0){next}

    ######  1  check Na+ adduct ####
    mz_index <- which(abs(peaks$mz[rt_index] - t_mz - NaAdduct) <= 0.01)

    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] - t_mz - NaAdduct))
      final_index_Na <- rt_index[mz_index]
      mz_Na <- peaks$mz[final_index_Na]
      Na_eic_matrix <- EIC_matrix(rawlcms, mz_Na, t_rt/60)
      Na_int_v <- peak_smooth(Na_eic_matrix[,2])

      Na_pearson <- round(cor(t_int_v, Na_int_v, method = "pearson"),2)
      if(Na_pearson >= 0.8) {
        if(peaks$Adduct[final_index_Na] != 0){
          peaks$Adduct[final_index_Na] <- paste0(peaks$Adduct[final_index_Na],
                                                 ";[M+Na]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_Na] <- paste0("[M+Na]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }
    ############ 2 check K+ adduct ###########
    mz_index <- which(abs(peaks$mz[rt_index] - t_mz - KAdduct) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index]- t_mz  - KAdduct))
      final_index_K <- rt_index[mz_index]
      mz_K <- peaks$mz[final_index_K]
      K_eic_matrix <- EIC_matrix(rawlcms, mz_K, t_rt/60)
      K_int_v <- peak_smooth(K_eic_matrix[,2])

      K_pearson <- round(cor(t_int_v, K_int_v, method = "pearson"),2)
      if(K_pearson >= 0.8) {
        if(peaks$Adduct[final_index_K] != 0){
          peaks$Adduct[final_index_K] <- paste0(peaks$Adduct[final_index_K],
                                                ";[M+K]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_K] <- paste0("[M+K]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }

    ############ 3 check NH4+ adduct ##########
    mz_index <- which(abs(peaks$mz[rt_index] - t_mz - NH4Adduct) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] - t_mz- NH4Adduct))
      final_index_NH4 <- rt_index[mz_index]
      mz_NH4 <- peaks$mz[final_index_NH4]
      NH4_eic_matrix <- EIC_matrix(rawlcms, mz_NH4 ,t_rt/60)
      NH4_int_v <- peak_smooth(NH4_eic_matrix[,2])

      NH4_pearson <- round(cor(t_int_v, NH4_int_v, method = "pearson"),2)
      if(NH4_pearson >= 0.8) {
        if(peaks$Adduct[final_index_NH4] != 0){
          peaks$Adduct[final_index_NH4] <- paste0(peaks$Adduct[final_index_NH4],
                                                  ";[M+NH4]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_NH4] <- paste0("[M+NH4]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }
    ############ 4.1 check loss H2O ######
    mz_index <- which(abs(peaks$mz[rt_index] + H2OAdduct - t_mz) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] + H2OAdduct - t_mz))
      final_index_H2O <- rt_index[mz_index]
      mz_H2O <- peaks$mz[final_index_H2O]
      H2O_eic_matrix <- EIC_matrix(rawlcms, mz_H2O,t_rt/60)
      H2O_int_v <- peak_smooth(H2O_eic_matrix[,2])
      H2O_pearson <- round(cor(t_int_v, H2O_int_v, method = "pearson"),2)
      if(H2O_pearson >= 0.8) {
        if(peaks$Adduct[final_index_H2O] != 0){
          peaks$Adduct[final_index_H2O] <- paste0(peaks$Adduct[final_index_H2O],
                                                  ";[M+H-H2O]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_H2O] <- paste0("[M+H-H2O]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }
    ############ 4.1 check addition H2O ######
    mz_index <- which(abs(peaks$mz[rt_index] - H2OAdduct - t_mz) <= 0.01)

    if(length(mz_index) !=0){
      #num_index <- which(peaks$cl[rt_index] == t_num)
      #rt_index <- rt_index[num_index]
      #if(length(rt_index)==0){next}

      mz_index <- which.min(abs(peaks$mz[rt_index] - H2OAdduct - t_mz))
      final_index_H2O <- rt_index[mz_index]
      mz_H2O <- peaks$mz[final_index_H2O]
      #mzRange_H2O <- c(mz_H2O-0.01, mz_H2O+0.01)
      # extract EIC for adduct, RT range: 40 seconds
      H2O_eic_matrix <- EIC_matrix(rawlcms, mz_H2O,t_rt/60)
      H2O_int_v <- peak_smooth(H2O_eic_matrix[,2])
      H2O_pearson <- round(cor(t_int_v, H2O_int_v, method = "pearson"),2)
      if(H2O_pearson >= 0.8) {
        #peaks$Adduct[i] <- paste0(peaks$Adduct[i],"[M+H]+_of_",
        #                                 peaks$ID_light[final_index_H2O]," ")
        if(peaks$Adduct[final_index_H2O] != 0){
          peaks$Adduct[final_index_H2O] <- paste0(peaks$Adduct[final_index_H2O],
                                                  ";[M+H+H2O]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_H2O] <- paste0("[M+H+H2O]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }

    ########### 5 check CH3OH adduct ##########
    # we believe the CH3OH adduct should be lower than the precursor #
    mz_index <- which(abs(peaks$mz[rt_index] - CH3OHAdduct - t_mz) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] - CH3OHAdduct - t_mz))
      final_index_CH3OH <- rt_index[mz_index]
      mz_CH3OH <- peaks$mz[final_index_CH3OH]
      CH3OH_eic_matrix <- EIC_matrix(rawlcms, mz_CH3OH,t_rt/60)
      CH3OH_int_v <- peak_smooth(CH3OH_eic_matrix[,2])
      CH3OH_pearson <- round(cor(t_int_v, CH3OH_int_v, method = "pearson"),2)
      if(CH3OH_pearson >= 0.8) {
        if(peaks$Adduct[final_index_CH3OH] != 0){
          peaks$Adduct[final_index_CH3OH] <- paste0(peaks$Adduct[final_index_CH3OH],
                                                    ";[M+H+CH3OH]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_CH3OH] <- paste0("[M+H+CH3OH]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }

    ########### 6 check [2M+H]+ dimer adduct #########
    # we believe the [2M+H]+ dimer adduct should be lower than the precursor #
    mz_index <- which(abs(peaks$mz[rt_index] - ( 2*t_mz - 1.007276)) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] - ( 2*t_mz - 1.007276)))
      final_index_dimer <- rt_index[mz_index]
      mz_dimer <- peaks$mz[final_index_dimer]
      dimer_eic_matrix <- EIC_matrix(rawlcms, mz_dimer,t_rt/60)
      dimer_int_v <- peak_smooth(dimer_eic_matrix[,2])
      dimer_pearson <- round(cor(t_int_v, dimer_int_v, method = "pearson"),2)
      if(dimer_pearson >= 0.8) {
        if(peaks$Adduct[final_index_dimer] != 0){
          peaks$Adduct[final_index_dimer] <- paste0(peaks$Adduct[final_index_dimer],
                                                    ";[2M+H]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_dimer] <- paste0("[2M+H]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }
    ########### 7 check [2M+NH4]+ dimer adduct #########
    mz_index <- which(abs(peaks$mz[rt_index] - ( 2*t_mz - 1.007276 + 17.026549)) <= 0.01)
    if(length(mz_index) !=0){
      mz_index <- which.min(abs(peaks$mz[rt_index] - ( 2*t_mz - 1.007276 + 17.026549)))
      final_index_dimerNH4 <- rt_index[mz_index]
      mz_dimerNH4 <- peaks$mz[final_index_dimerNH4]
      dimerNH4_eic_matrix <- EIC_matrix(rawlcms, mz_dimerNH4,t_rt/60)
      dimerNH4_int_v <- peak_smooth(dimerNH4_eic_matrix[,2])
      dimerNH4_pearson <- round(cor(t_int_v, dimerNH4_int_v, method = "pearson"),2)
      if(dimerNH4_pearson >= 0.8) {
        if(peaks$Adduct[final_index_dimerNH4] != 0){
          peaks$Adduct[final_index_dimerNH4] <- paste0(peaks$Adduct[final_index_dimerNH4],
                                                       ";[2M+NH4]+_of_S",samNum,"_",peaks$featureID[k])
        }else{peaks$Adduct[final_index_dimerNH4] <- paste0("[2M+NH4]+_of_S",samNum,"_",peaks$featureID[k])}
      }
    }
  }
  adduct <- peaks[peaks$Adduct != 0,]
  colnames(peaks) <- replace(colnames(peaks),colnames(peaks)== "maxo", "Intensity")
  return(peaks)
}


#' Title: identify isotopes #
#'
#' @param chlorine_tb data frame, Cl peaks summarized in a table
#' @param rawlcms xcms object processed by xcmsRaw
#' @param samNum numeric
#' @param rt_tol numeric, retention time in seconds
#' @param mz_tol numeric, m/z  in Da
#'
#' @return data frame with isotope labelled
#' @import xcms
#'

findIsotope <- function(chlorine_tb, rawlcms, samNum =1 , rt_tol =10, mz_tol=0.01){
  peaks <- chlorine_tb
  peaks$isotope <- 0

  for(k in 1:nrow(peaks)){
    t_mz <- peaks$mz[k]
    t_rt <- peaks$rt[k]
    t_num <- peaks$cl[k]
    t_eic_matrix <- EIC_matrix(rawlcms, t_mz, t_rt/60)
    t_int_v <- t_eic_matrix[,2]
    t_int_v <- peak_smooth(t_int_v)

    rt_index <- which(abs(peaks$rt - t_rt) <= rt_tol)
    if(length(rt_index) == 0){next}

    # detect 13C
    C13_mz_index <- which(abs(peaks$mz[rt_index] - t_mz - 1.003355) <= 0.01)
    if(length(C13_mz_index) != 0){
      if(length(C13_mz_index) >1){
        C13_mz_index <- which.min(abs(peaks$mz[rt_index] - t_mz- 1.003355))
      }
      final_c13_index <- rt_index[C13_mz_index]

      c13_mz <- peaks$mz[final_c13_index]
      c13_eic <- EIC_matrix(rawlcms,c13_mz,t_rt/60)
      c13_eic[,2] <- peak_smooth(c13_eic[,2])
      cor <- cor(t_int_v, c13_eic[,2])
      if(cor >= 0.8){
        if(peaks$isotope[final_c13_index] != 0){
          peaks$isotope[final_c13_index] <- paste0(peaks$isotope[final_c13_index],
                                                   ";M+1_for_S",samNum,"_",peaks$featureID[k],
                                                   collapse = ";")
        }else{
          peaks$isotope[final_c13_index] <- paste0("M+1_for_S",samNum,"_",peaks$featureID[k])
        }
      }
    }

    # detect 37Cl
    Cl37_mz_index <- which(abs(peaks$mz[rt_index] - t_mz - 1.99705) <= 0.01)
    if(length(Cl37_mz_index) != 0){
      if(length(Cl37_mz_index) >1){
        Cl37_mz_index <- which.min(abs(peaks$mz[rt_index] - t_mz-1.99705))
      }
      final_cl37_index <- rt_index[Cl37_mz_index]

      cl37_mz <- peaks$mz[final_cl37_index]
      cl37_eic <- EIC_matrix(rawlcms,cl37_mz,t_rt/60)
      cl37_eic[,2] <- peak_smooth(cl37_eic[,2])
      cor <- cor(t_int_v, cl37_eic[,2])
      if(cor >= 0.8){
        if(peaks$isotope[final_cl37_index] != 0){
          peaks$isotope[final_cl37_index] <- paste0(peaks$isotope[final_cl37_index],
                                                    ";M+2_for_S",samNum,"_",peaks$featureID[k],
                                                    collapse = ";")
        }else{
          peaks$isotope[final_cl37_index] <- paste0("M+2_for_S",samNum,"_",peaks$featureID[k])
        }
      }
    }
  }
  return(peaks)
}

#' Title: find in source fragments
#'
#' @param peaks data frame, all peaks
#' @param chlorine_tb data frame, all Cl-containing peaks
#' @param rawfile_dir string, file folder storing the raw mzML/mzXML files
#' @param lcmsfile string, mzML/mzXML file
#' @param sampleNumber numeric, sample number
#' @param ISFtable data frame, in source fragment from all peaks
#'
#' @return data frame with in source fragment labelled
#' @export
#' @import xcms

findInsourcefrag <- function(peaks, chlorine_tb, rawfile_dir, lcmsfile, sampleNumber = 1,
                             ISFtable){
  result <- ISFtable
  isf_index <- which(result$ISF_level !=0)
  if(length(isf_index)  == 0){
    chlorine_tb$ISF_level <- 0
    return(chlorine_tb)
  }

  for(m in isf_index){
    # if(result$Num_Level2[m] == 0 & result$Num_Level1[m] == 0)
    # only keep the in source fragment level 1
    levels <- strsplit(result$ISF_level[m], ";")[[1]]
    levels <- unique(levels)
    levels <- levels[grepl("Level_1", levels)]

    if(length(levels) == 0){
      result$ISF_level[m] <- 0
    }else{
      precursor_id <- unlist(strsplit(levels, "<-Level_1"))
      cl_precursor_id <- precursor_id[which( precursor_id %in% result$featureID)]

      if(length(cl_precursor_id) ==0){
        # all precursors are not cl, the ISFrag is not good
        result$ISF_level[m] <- 0
      }else{
        result$ISF_level[m] <- paste0(cl_precursor_id, collapse = ";")
      }
    }

  }

  chlorine_tb$ISF_level <- result$ISF_level
  potential_isfrag_i <- which(chlorine_tb$ISF_level != 0)

  if(length(potential_isfrag_i) == 0){
    return(chlorine_tb)
  }

  true_pos_tb <- chlorine_tb[chlorine_tb$Adduct ==0,]
  true_pos_tb <- true_pos_tb[true_pos_tb$isotope ==0,]
  true_id <- true_pos_tb$featureID

  for(k in 1:length(potential_isfrag_i)){
    i <- potential_isfrag_i[k]
    precursor_i <- strsplit(chlorine_tb$ISF_level[i],";")[[1]]
    final_precursor <- precursor_i[which(precursor_i %in% true_id)]

    # if the precursor is not the adduct form or isotope form
    # the in source fragment identification could be true
    if(length(final_precursor) == 0){
      chlorine_tb$ISF_level[i] <- 0
    }else{
      #break
      final_precursor <- paste0( "S",sampleNumber,"_", final_precursor)
      chlorine_tb$ISF_level[i] <- paste0( final_precursor ,"<-Level_1", collapse = ";")
    }
  }
  return(chlorine_tb)
}
