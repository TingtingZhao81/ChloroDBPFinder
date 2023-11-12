#' Title
#'
#' @param featureTable data frame, peaks
#' @param Cl_db database, check "db_demo"
#' @param ion_mode string, ionization mode "P" for positive mode
#' @param ref_mz_tol numeric, mass tolerance in Da
#' @param dp_score MS2 spectra similarity score threshold
#' @param dp_num number of matched fragments
#'
#' @return data frame
#' @export
#'
annonateFeature <- function( featureTable, Cl_db, ion_mode = "P", ref_mz_tol =25,
                          dp_score = 70, dp_num = 2){
  featureTable$name <- 0
  featureTable$match_num <- 0
  featureTable$score <- 0

  Cl_db$Precursor_mz <- as.numeric(Cl_db$Precursor_mz )
  table(Cl_db$Ion_mode)
  Cl_db <- Cl_db[Cl_db$Ion_mode == ion_mode,]

  annotation_num <- 0
  for(i in 1:nrow(featureTable)){
    if(featureTable$ms2_mz[i] == 0 | featureTable$ms2_mz[i] == "" | is.na(featureTable$ms2_mz[i] )){next}
    t_mz <- featureTable$mz[i]
    index <- which(abs(Cl_db$Precursor_mz - t_mz) <= t_mz*ref_mz_tol * 0.000001)
    ms2_mz  <- as.numeric(strsplit(featureTable$ms2_mz[i],";")[[1]])
    ms2_int <- as.numeric(strsplit(featureTable$ms2_int[i],";")[[1]])
    if(length(index) != 0){
      score_v <- c()
      match_num <- c()
      for(k in index){
        db_ms2mz <- as.numeric(strsplit(Cl_db$MS2mz[k],";")[[1]])
        db_ms2int <- as.numeric(strsplit(Cl_db$MS2int[k],";")[[1]])
        dp <- dotProduct(featureTable$ms2_mz[i],featureTable$ms2_int[i],Cl_db$MS2mz[k], Cl_db$MS2int[k])
        score_v <- c(score_v, as.numeric(dp[1]))
        match_num <- c(match_num, as.numeric(dp[2]))
      }
      idx <- intersect(which(score_v >= dp_score) ,which(match_num >= dp_num))
      if(length(idx) == 0 ){
        next
      }
      annotation_num <- annotation_num + 1

      order <- order(score_v[idx], decreasing = TRUE)

      featureTable$name[i] <- paste0(Cl_db$Name[index[idx]][order], collapse = ";")
      featureTable$score[i] <- paste0(round(score_v[idx][order])/100, collapse = ";")
      featureTable$match_num[i] <- paste0(match_num[idx][order], collapse = ";")
    }
  }

  return(featureTable)
}
