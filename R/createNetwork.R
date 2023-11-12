
#' Title: establish network
#'
#' @param featureTable data frame, all peaks
#' @param seed_tb data frame, reaction precursors or confidently annotated features
#' @param reaction_pathways data frame, possible reaction pathways
#' @param reaction_pathway_mz_tol numeric, m/z tolerance in reaction pathway
#' @param nw_spectra_score numeric, GNPS spectra similarity threshold
#' @param nw_spectra_match_num numeric, least number of matched fragments
#'
#' @return list, all network and explainable network
#' @export
#'
createNetwork <- function(featureTable, seed_tb, reaction_pathways, reaction_pathway_mz_tol = 0.01,
                          nw_spectra_score=0.5, nw_spectra_match_num=3) {
  target_l <- seed_tb
  ft <- featureTable
  # add extra column into cl_tb
  append_col <- colnames(target_l)[which( !colnames(target_l)  %in% colnames(ft))]

  for(i in 1:length(append_col)){
    ft <- cbind(ft,"")
    colnames(ft)[ncol(ft)] <- append_col[i]
  }

  #### merge the target_l  and identified cl-tb
  index <- which( !ft$featureID %in% target_l$featureID)
  ft <- ft[index,]
  ft$label <- "unknown"
  merged_tb <- target_l
  colnames(target_l)
  col_indexes <- c()
  for(i in colnames(target_l)){
    col_indexes <- append(col_indexes, which(colnames(ft) == i ), length(col_indexes) )
  }
  merged_tb <- rbind(merged_tb, ft[,col_indexes])
  merged_tb$cl <- as.numeric(merged_tb$cl)

  colnames_nw <- c("node1", "label1","mz1", "rt1","cl1", "annotation1", "mf1", "ms2_mz1", "ms2_int1",
                   "node2", "label2","mz2", "rt2","cl2", "annotation2", "mf2", "ms2_mz2", "ms2_int2",
                   "mzdiff", "reaction_pathway","change_of_mf" ,"spectral_score","spectral_num")
  network_tb <- as.data.frame(matrix(ncol=length(colnames_nw), nrow=nrow(merged_tb)*(nrow(merged_tb) -1)/2))
  colnames(network_tb) <- colnames_nw

  track_row <- 1
  cols_infor <- vector()
  names <- c("featureID","label", "mz", "rt", "cl", "annotation", "MF", "ms2_mz", "ms2_int")
  for(l in 1:length(names)){
    cols_infor <- append(cols_infor, which(colnames(merged_tb) == names[l]), length(names) )
  }

  for(i in 1:(nrow(merged_tb) -1)){

    ms2_mz1 <- merged_tb$ms2_mz[i]
    ms2_int1 <- merged_tb$ms2_int[i]

    start_i <- i+1
    for(j in start_i:nrow(merged_tb)){
      network_tb[track_row, 1:length(names)] <- merged_tb[i, cols_infor]
      network_tb[track_row, (length(names)+1):(length(names)*2)] <- merged_tb[j, cols_infor]

      ### check reaction pathway ###
      mzdiff <- merged_tb$mz[j] - merged_tb$mz[i]
      network_tb$mzdiff[track_row] <- mzdiff

      base_dif <- ( merged_tb$cl[j] - merged_tb$cl[i]) * 33.96102
      mz_diff <- mzdiff - base_dif

      index <- which(abs(reaction_pathways$mass_change - mz_diff)  <= reaction_pathway_mz_tol)
      if(length(index) != 0){
        # the mass difference match certain reaction pathway
        path <- paste0(reaction_pathways$reaction_pathway[index], collapse = ";")
        network_tb$reaction_pathway[track_row] <- paste0(merged_tb$cl[j] - merged_tb$cl[i], "Cl_", path)
        mf_change <- paste0(reaction_pathways$formula_change[index], collapse = ";")
        network_tb$change_of_mf[track_row] <- paste0(merged_tb$cl[j] - merged_tb$cl[i], "Cl_", mf_change)
      }



      ### check spectral similarity ###
      ms2_mz2 <- merged_tb$ms2_mz[j]
      ms2_int2 <- merged_tb$ms2_int[j]

      if(ms2_mz1 != "" & ms2_mz2 != ""){
        ms2_1 <- as.data.frame(cbind(as.numeric(strsplit(ms2_mz1, split = ";")[[1]]),
                                     as.numeric(strsplit(ms2_int1, split = ";")[[1]])))
        ms2_2 <- as.data.frame(cbind(as.numeric(strsplit(ms2_mz2, split = ";")[[1]]),
                                     as.numeric(strsplit(ms2_int2, split = ";")[[1]])))

        gnps <- GNPS_score(ms2_1, merged_tb$mz[i] ,ms2_2, merged_tb$mz[j] )

        network_tb$spectral_score[track_row] <- gnps$score
        network_tb$spectral_num[track_row] <- gnps$num
      }

      track_row <- track_row + 1
    }
  }

  similarity_nw <- network_tb[!is.na(network_tb$spectral_score),]
  similarity_nw <- similarity_nw[similarity_nw$spectral_score >= nw_spectra_score,]
  similarity_nw <- similarity_nw[similarity_nw$spectral_num >= nw_spectra_match_num,]
  explainable_nw <- similarity_nw[!is.na(similarity_nw$reaction_pathway),]

  return(list(network_tb, explainable_nw))
}



#' Title: dot product function
#' exp_fragMz, exp_fragInt, ref_fragMz, ref_fragMzInt: separate by ";"
#'
#' @param exp_fragMz string, experimental fragment m/z separated by ";"
#' @param exp_fragInt string, experimental fragment intensity separated by ";"
#' @param ref_fragMz string,  reference fragment m/z separated by ";"
#' @param ref_fragMzInt string,  reference fragment intensity separated by ";"
#' @param ms2_tol numeric, m/z tolerance in Da
#'
#' @return list, dp score, count of matched fragment, mz
#'
dotProduct <- function(exp_fragMz, exp_fragInt, ref_fragMz, ref_fragMzInt, ms2_tol =0.02) {
  # get the vector of experimental fragMz and fragInt
  exp_fragMz <- as.numeric(strsplit(exp_fragMz, ";")[[1]])
  exp_fragInt <- as.numeric(strsplit(exp_fragInt, ";")[[1]])
  exp_fragInt <- 100*exp_fragInt/max(exp_fragInt)
  # calculate the length of intensity vector
  A <- sqrt(sum(exp_fragInt^2))
  # record how many experimental fragment
  totalFrag <- length(exp_fragMz)
  # find the reference frag_mz and frag_int
  ref_fragMz <- as.numeric(strsplit(ref_fragMz, ";")[[1]])
  ref_fragInt <- as.numeric(strsplit(ref_fragMzInt, ";")[[1]])
  ref_fragInt <- 100*ref_fragInt/max(ref_fragInt)
  # calculate the length of intensity vector
  B <- sqrt(sum(ref_fragInt^2))
  # create the aligned data frame
  int_dataframe <- data.frame(matrix(nrow = length(ref_fragMz),ncol=4))
  colnames(int_dataframe) <- c("ref_fragMz","ref_fragInt", "exp_fragInt","exp_fragMz")
  int_dataframe[,1] <- as.data.frame(ref_fragMz)
  int_dataframe[,2] <- as.data.frame(ref_fragInt)
  count <- 0 # record how many counts found
  for (j in 1:length(ref_fragMz)){
    fragMz_index <- which( exp_fragMz %in% exp_fragMz[abs(exp_fragMz - int_dataframe[j,1]) <= ms2_tol] )
    # no frag found in t_fragment, assign 0
    if ( length(fragMz_index)==0) {
      int_dataframe[j,3] <- 0
      next
    }
    # assign the found intensity to dataframe
    frag_index <- fragMz_index[which.max(exp_fragInt[fragMz_index])]
    int_dataframe[j,3] <- exp_fragInt[frag_index[1]]
    int_dataframe[j,4] <- exp_fragMz[frag_index[1]]
    count <- count + 1
    #remove the matched fragment from experimental
    exp_fragMz <- exp_fragMz[-frag_index]
    exp_fragInt <- exp_fragInt[-frag_index]
  }
  # calculate the dp
  dp <- sum(int_dataframe$ref_fragInt * int_dataframe$exp_fragInt)
  mz <- paste0(int_dataframe$ref_fragMz[complete.cases(int_dataframe$exp_fragMz)],
               collapse = ";")
  score <- 100*dp/(A*B)
  results <- c(score, count, mz)
  # return the number of matched fragments and dot product score
  return(results)
}


#' Title: GNPS spectral similarity function #####
#'  library: clue
#'  input: x,y: data frame containing m/z and intensity in ms2 spectra;
#          x.pre, y.pre: precursor ion m/z'
#' @param x data frame, ms2 spectra with m/z column and intensity column
#' @param x.pre numeric, precursor ion m/z
#' @param y data frame, ms2 spectra with m/z column and intensity column
#' @param y.pre numeric, precursor ion m/z
#' @param mz.tol numeric, m/z tolerance
#'
#' @return list, score and matched number
#' @import clue
#' @importFrom stats complete.cases
#'
GNPS_score <- function(x, x.pre, y, y.pre, mz.tol=0.02){
  colnames(x) <- c("mz","int")
  colnames(y) <- c("mz","int")
  #delete pre
  score <- 0
  num <- 0
  remove_index <- which( abs(x$mz -x.pre) <= 0.02)
  if (length(remove_index) !=0)  {
    x <- x[-remove_index,]
    if(nrow(x) == 0) {return(as.data.frame(cbind(score,num)))}
  }

  remove_index <- which( abs(y$mz -y.pre) <= 0.02)
  if (length(remove_index) !=0) {
    y <- y[-remove_index,]
    if(nrow(y) == 0) {return(as.data.frame(cbind(score,num)))}
  }

  dm <- y.pre - x.pre
  #square root transforms
  x[,2] <- sqrt(x[,2])
  y[,2] <- sqrt(y[,2])
  #scale intensity values
  x[,2] <- x[,2]/sqrt(sum(x[,2]^2))
  y[,2] <- y[,2]/sqrt(sum(y[,2]^2))
  # create alignment matrix
  alignment <- data.frame(matrix(ncol = 5))
  colnames(alignment) <- c("mz.x", "int.x", "int.y","mz.y","note")
  h <- 1
  for (m in 1:nrow(x)){
    mz.diff1 <- abs(x[m,1] - y[,1])
    mz.diff2 <- abs(x[m,1] - y[,1] + dm)
    if(min(mz.diff1) <= mz.tol){
      alignment[h,1] <- x[m,1]
      alignment[h,2] <- x[m,2]
      alignment[h,3] <- y[mz.diff1 == min(mz.diff1),2][1]
      alignment[h,4] <- y[mz.diff1 == min(mz.diff1),1][1]
      alignment[h,5] <- "direct"
      h <- h + 1
    }
    if(dm != 0){
      if(min(mz.diff2) <= mz.tol){
        alignment[h,1] <- x[m,1]
        alignment[h,2] <- x[m,2]
        alignment[h,3] <- y[mz.diff2 == min(mz.diff2),2][1]
        alignment[h,4] <- y[mz.diff2 == min(mz.diff2),1][1]
        alignment[h,5] <- "neutralLoss"
        h <- h + 1
      }
    }
  }
  alignment <- alignment[complete.cases(alignment),]

  # create a matrix for selecting the highest scoring subset of matching peaks using 'clue' library
  # each peak is matched to at most one peak in another spectrum
  if(nrow(alignment)==0){
    score <- 0
    num <- 0
  }
  if(nrow(alignment)>0){
    mzfrag.x <- unique(alignment$mz.x)
    mzfrag.y <- unique(alignment$mz.y)
    max.length <- max(length(mzfrag.x),length(mzfrag.y))
    num <- min(length(mzfrag.x),length(mzfrag.y))
    matrix <- matrix(0, nrow = max.length, ncol = max.length) #matrix nrow=ncol
    #name the matrix rows and cols with mz
    if(length(mzfrag.x)<=length(mzfrag.y)){
      rownames(matrix) <- c(mzfrag.x, rep(0,(max.length-length(mzfrag.x))))
      colnames(matrix) <- mzfrag.y
    }
    if(length(mzfrag.x)>length(mzfrag.y)){
      rownames(matrix) <- mzfrag.x
      colnames(matrix) <- c(mzfrag.y,rep(0,(max.length-length(mzfrag.y))))
    }

    #fill the matrix with intensity product
    for(h in 1:nrow(alignment)){
      matrix[mzfrag.x==alignment[h,1],mzfrag.y==alignment[h,4]] <- alignment[h,2] * alignment[h,3]
    }
    if(length(mzfrag.x)<length(mzfrag.y)){matrix[(length(mzfrag.x)+1):max.length,]<-0}
    if(length(mzfrag.x)>length(mzfrag.y)){matrix[,(length(mzfrag.y)+1):max.length]<-0}
    #LSAP problem in 'clue' library
    optimal <- clue::solve_LSAP(matrix, maximum = TRUE)
    #calculate GNPS score
    GNPS <- sum(matrix[cbind(seq_along(optimal),optimal)])
    score <- as.numeric(GNPS)
  }
  t <- as.data.frame(cbind(score, num))
  return(t)
}


