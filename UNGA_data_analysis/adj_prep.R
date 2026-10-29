library(parallel)
library(dplyr)
library(tidyr)
library(RSpectra)
library(cluster)
library(Rfast)
library(data.table)
library(readr)
library(irlba)
library(softImpute)
library(Matrix)
library(pROC)
library(IMIFA)
library(rlist)
library(rdist)
library(latentnet)
library(network)
library(unvotes)
library(ergm)

##### similarity functions: cosine and jaccard ##################

#cosine similarity #
normalize <- function(X){
  ## sometimes it may happen some countries give NV to all resolutions.
  ## in that denom rowsums(X^2) for that row becomes 0, hence it produces NaN.
  ## For example, JAPAN in 1956 gives NV to all resolutions.
  ## so we put 0 to all.
  row_sums_sq <- rowSums(X^2)  
  non_zero_rows <- which(row_sums_sq > 0)  
  X[non_zero_rows, ] <- X[non_zero_rows, ] / sqrt(row_sums_sq[non_zero_rows])  # Normalize only nonzero rows
  return(X)
}

# jaccard similarity #
jaccard_similarity <- function(vector_1,vector_2){
  if(length(vector_1) == 0 || length(vector_2) == 0){
    
    stop("error! one of the vectors can not be empty")
    
  }else if (length(vector_1) != length(vector_2)) {
    
    stop("error! lengths of the vectors do not match.")
    
  }else{
    
    matches <- vector_1 == vector_2
    intersection <- vector_1[matches]
    jaccard_sim = length(intersection)/(length(vector_1)+length(vector_2)-length(intersection))
    
  }
  return(jaccard_sim)
}



############# adj prep from package unvotes ####################

adj_prep_from_unga_r <- function(input_year, similarity_measure, threshold = c('q1','q2','q3','mean')){
  
  un_roll_calls = data.frame(unvotes::un_roll_calls)
  un_roll_call_issues = data.frame(unvotes::un_roll_call_issues)
  un_votes = data.frame(unvotes::un_votes)
  
  d1 = full_join(un_roll_calls, un_roll_call_issues, by = "rcid")
  d2 = inner_join(d1, un_votes, by = "rcid",relationship = "many-to-many")
  
  d2_selected = d2[year(d2$date) == input_year,]
  
  #names(d2_selected)[names(d2_selected) == "country_code"] <- "iso2c"
  #unique(d2_selected$vote) ## yes,abstain, no
  df = d2_selected

  vote_mapping <- c('yes' = 1, 'no' = -1, 'abstain' = 0)
  
  df$vote_numeric <- vote_mapping[df$vote]
  
  df_final <-  df[c('rcid','descr','country','vote_numeric')]
  
  ## important as 1997 description contain ""
  df_final <- df_final[df_final$descr != "", ]
  
  countries <- sort(unique(df_final$country))
  descriptions <- sort(unique(df_final$descr))

  # Create an empty matrix filled with 0
  vote_matrix <- matrix(0, 
                        nrow = length(countries),
                        ncol = length(descriptions),
                        dimnames = list(countries, descriptions))
  
  # Fill matrix
  for (i in seq_len(nrow(df_final))) {
    row_country <- df_final$country[i]
    col_descr   <- df_final$descr[i]
    value       <- df_final$vote_numeric[i]
    
    vote_matrix[row_country, col_descr] <- value
  }
  
  rownames(vote_matrix) <- countries
  
  #####################################################
  ############ data preparation ends ##################
  ################ ************ #######################
  
  ###################################################################################
  ############ prepare the adj matrix ###############################################
  ############ construct adj using similarity between two rows of vote_matrix #######
  #### construct binary adj from wtd adj using different thresholds: preference-q2 ##
  ##################################################################### #############
  
    
  if (similarity_measure == 'cosine'){
    
    normalized_matrix = normalize(vote_matrix)
    
    Adj <- round(as.matrix(tcrossprod(normalized_matrix)),3)# without rounding it may violate symmetricity due to R aprroximation
    
  }else if (similarity_measure == 'jaccard'){
    
    Adj = outer(1:nrow(vote_matrix), 1:nrow(vote_matrix), Vectorize(function(i, j) jaccard_similarity(vote_matrix[i, ], vote_matrix[j, ])))
    
    Adj = round(as.matrix(Adj),3)
  }else {
    
    stop("Error!choose cosine or jaccard similarity")
  }
  
  ## no self-loop
  diag(Adj)<- 0
  Adj[Adj<0]<-0  # making the negative weights = 0, because cosine similarity \in [-1,1]
  
  ## making it binary based on threshold
  if (threshold == 'q1'){
    Adj_binary = (Adj > quantile(Adj,0.25)) * 1
  }else if (threshold == 'q2'){
    Adj_binary<- (Adj > median(Adj)) * 1 #### robust choice of threshold ####
  }else if (threshold == 'q3'){
    Adj_binary = (Adj > quantile(Adj,0.75)) * 1
  }else{
    Adj_binary = (Adj > mean(Adj)) * 1
  }
  
  return(Adj_binary)
  
}

