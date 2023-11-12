#' Title: align the features across samples
#'
#' @param file_dir string, the folder storing the feature table to align
#' @param filePattern string, common file names of the feature tables
#' @param align_mz_tol numeric, m/z tolerance in Dalton
#' @param align_rt_tol numeric, rt tolerance in seconds
#' @param align_num logic, true or false when align the features
#'
#' @importFrom utils read.csv
#' @return data frame
#' @export
#'
alignFeature <- function(file_dir, filePattern, align_mz_tol = 0.01, align_rt_tol = 20, align_num =FALSE){
  cl_tbs <- list.files(pattern = filePattern, path = file_dir)
  if(length(cl_tbs)!= 1){
    # do alignment
    colname <- c("mz","rt","cl", "IDs","ms2_mz","ms2_int", "ms2_origin","freq") # IDs are recorded in the fourth column
    colnames <- colname
    for(k in 1:length(cl_tbs)){
      sampleName <- strsplit(cl_tbs[k],split = filePattern)[[1]][1]
      colnames <- c(colnames,sampleName)
    }
    align_tb <- as.data.frame(matrix(ncol=length(colnames)))
    colnames(align_tb) <- colnames

    tb1 <- read.csv(paste0(file_dir ,"/",cl_tbs[1]))
    #tb1 <- tb1[tb1$cor >=cor_thresh,]
    message(nrow(tb1), "_cl_in_", strsplit(cl_tbs[1],split = filePattern)[[1]][1])
    colnames(tb1)
    indi <- align_tb
    for(k in 1:nrow(tb1)){
      indi$mz <- tb1$mz[k]
      indi$rt <- tb1$rt[k]
      indi$cl <- tb1$cl[k]
      indi$IDs <- tb1$featureID[k]
      indi$ms2_mz <- tb1$ms2_mz[k]
      indi$ms2_int <- tb1$ms2_int[k]
      if(!(tb1$ms2_mz[k] == 0 | tb1$ms2_mz[k] == "" | is.na(tb1$ms2_mz[k] ))){
        indi$ms2_origin <- 1
      }
      indi[,length(colname)+1] <- tb1$Intensity[k]
      align_tb <- rbind(align_tb, indi)
    }
    align_tb <- align_tb[-1,]

    for(k in 2:length(cl_tbs)){
      # k th files
      tb <- read.csv(paste0(file_dir ,"/",cl_tbs[k]))
      #tb <- tb[tb$cor >= cor_thresh,]
      message(nrow(tb), "_cl_in_",
              strsplit(cl_tbs[k],split = "_")[[1]][1])
      # align with exiting cl-compounds
      for(j in 1:nrow(tb)){
        t_mz <- tb$mz[j]
        t_rt <- tb$rt[j]
        t_num <- tb$cl[j]
        if(align_num){
          # do consider the Cl number match #
          num_index <- which(align_tb$cl == t_num)
        }else{
          # do not consider the Cl number match
          num_index <- 1:nrow(align_tb)
        }

        rt_index <- which(abs(align_tb$rt[num_index] - t_rt)
                          <= align_rt_tol)
        mz_index <- which(abs(align_tb$mz[num_index][rt_index] - t_mz)
                          <= align_mz_tol)

        if(length(mz_index) == 0){
          # not found in the existing table, add a new row
          indi <- as.data.frame(matrix(ncol=length(colnames)))
          colnames(indi) <- colnames
          indi$mz <- t_mz
          indi$rt <- t_rt
          indi$cl <- t_num
          indi$IDs <- tb$featureID[j]
          indi$ms2_mz <- tb$ms2_mz[j]
          indi$ms2_int <- tb$ms2_int[j]
          if( !(tb$ms2_mz[j] == 0 | tb$ms2_mz[j] == "" | is.na(tb$ms2_mz[j] ))){
            indi$ms2_origin <- k
          }
          indi[,length(colname) +k] <- tb$Intensity[j]
          align_tb <- rbind(align_tb, indi)

        }else{
          # if match multiple features in the aligned table
          # evaluate base on the score
          if(length(mz_index) > 1){
            mz_score_v <- 1 - abs(align_tb$mz[rt_index[mz_index]] - t_mz)/align_mz_tol

            rt_score_v <- 1 - abs(align_tb$rt[rt_index[mz_index]] - t_rt)/align_rt_tol

            mz_index <- mz_index[which.max(mz_score_v + rt_score_v)]
          }

          # found in the existing table
          final_index <- num_index[rt_index[mz_index]]
          # update the ID
          align_tb$IDs[final_index] <- paste0(align_tb$IDs[final_index] ,";",tb$featureID[j])
          # update the number of Cl if do not consider the number of Cl in the alignment
          # take the majority to decide the number of Cl in the compound
          if(!align_num){
            align_tb$cl[final_index] <- paste0(align_tb$cl[final_index] ,";",tb$cl[j])
          }
          align_tb[final_index,length(colname)+k] <- tb$Intensity[j]
          # check MS2 whether need to update
          # has ms2
          # 2. check the intensity, larger than record one
          if( !(tb$ms2_mz[j] == 0 | tb$ms2_mz[j] == "" | is.na(tb$ms2_mz[j] ))){
            if(align_tb$ms2_mz[final_index] == 0 | align_tb$ms2_mz[final_index] == "" | is.na(align_tb$ms2_mz[final_index])){
              # no ms2 , add ms2 into the table
              align_tb$ms2_origin[final_index] <- k
              align_tb$ms2_mz[final_index] <- tb$ms2_mz[j]
              align_tb$ms2_int[final_index] <- tb$ms2_int[j]
            }else{
              # has ms2
              int <- align_tb[final_index, length(colname) + align_tb$ms2_origin[final_index]]
              if(int < tb$Intensity[j]){
                #update
                align_tb$ms2_origin[final_index] <- k
                align_tb$ms2_mz[final_index] <- tb$ms2_mz[j]
                align_tb$ms2_int[final_index] <- tb$ms2_int[j]
              }
            }
          }

        }
      }
    }
  }

  #write.csv(align_tb, "alignment before merge the predited Cl numbers.csv", row.names = FALSE)
  #align_tb <- read.csv("alignment before merge the predited Cl numbers.csv")
  start_n <-  ncol(align_tb) - length(cl_tbs) + 1
  end_n <- ncol(align_tb)

  for(i in 1:nrow(align_tb)){
    align_tb$freq[i] <- length(strsplit(align_tb$IDs[i], ";")[[1]])

    # more than two samples have and do not use number of cl to align features
    # then we need to use the majority vote to decide the number Cl
    if(align_tb$freq[i] != 1 && !align_num){

      # if not consider the number Cl in the alignment, then take the majority vote
      cl_v <- as.numeric(strsplit(align_tb$cl[i],";")[[1]])
      cl_tb <- as.data.frame(table(cl_v))
      # cl_tb$cl_v <- as.factor(cl_tb$cl_v)

      num_cl <- cl_tb$cl_v[which(cl_tb$Freq == max(cl_tb$Freq)) ]
      if(length(num_cl) == 1){
        align_tb$cl[i] <- as.character(num_cl)
      }else{
        # if number of votes are equal for different Cl number
        # use the number of cl from the highest intensity
        int_v <- vector()
        for(l in start_n:end_n){
          if(!is.na(align_tb[i, l]) || length(is.na(align_tb[i, l])) == 0){
            int_v <- c(int_v, align_tb[i, l])
          }
        }
        align_tb$cl[i] <- cl_v[which.max(int_v)]
      }
    }
  }
  return(align_tb)
}
