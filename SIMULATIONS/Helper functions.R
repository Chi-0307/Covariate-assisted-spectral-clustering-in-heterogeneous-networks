######### Helper functions ####

###=== Descriptions of the code =======###
# 1. Contains functions to be used inside hetcov,homcov,het and hom clustering.
# for ex - generate adj, calculating misclustering error rate etc. No of types has been hardcoded to T=2
# 2. Contains four main functions that computes no of misclustered nodes for both the types corres.to hetcov,homcov
# het and hom clusterings.
################################################################################

library(RSpectra)
library(clue)  # For the 'solve_LSAP' function to find optimal permutations


## generate adjacency ##
gen_adj<-function(membership_mat,block_prob_mat){
  #expected_A = membership_mat %*% block_Prob_mat %*% t(membership_mat) 
  expected_A = tcrossprod(membership_mat %*% block_prob_mat, membership_mat)
  N = nrow(membership_mat)
  #Generate the adjacency matrix A from expected_A (E(A|Z)) as parameter
  A = matrix(rbinom(N * N, 1, expected_A), nrow = N, ncol = N)
  #Set diagonal elements to 0 (no self-loops)
  diag(A) <- 0
  A[lower.tri(A)] = t(A)[lower.tri(A)] ### making adj symmetric
  return(A)
}

### generate regularized laplacian###
gen_regularized_laplacian <- function(A){
  degrees = rowSums(A)
  tau = mean(degrees+0.0001) # for stabilizing
  D_tau_power_neg_half = diag(1/sqrt(degrees+tau))
  L_tau <-  tcrossprod(crossprod(D_tau_power_neg_half, A), D_tau_power_neg_half)
  return(L_tau)
}


#get top eigen vector matrix of the Laplacian 
get_eigenvector_matrix <- function(L, num_clusters) {
  eigen_decomp <- eigen(L)
  # Select the eigenvectors corresponding to the largest(in abs value) nonzero eigenvalues
  eigenvalues <- eigen_decomp$values
  eigenvectors <- eigen_decomp$vectors
  
  #Sort eigenvalues by absolute value in descending order
  sorted_indices <- order(abs(eigenvalues), decreasing = TRUE)
  
  # Select the top T*K eigenvalues and corresponding eigenvectors
  top_indices <- sorted_indices[1:num_clusters]
  
  ## finally get the eigenvector matrix 
  V <- eigenvectors[, top_indices]
  return(V)
}


### not needed ####
#get within sum of sq for type_1,type_2 nodes separately
wss_alpha <- function(alpha,L_tau_sq,X,T,K,n) {
  
  L_tilde = L_tau_sq + alpha * tcrossprod(X)
  V <- RSpectra::eigs_sym(L_tilde, T*K, "LM")$vectors
  
  V_type1 = V[1:(n/2),]
  V_type2 = V[((n/2)+1):n,]
  
  # Apply k-means clustering to the eigenvectors
  k_means_output_type1 <- kmeans(V_type1, centers = K, iter.max = 100, nstart = 20)
  k_means_output_type2 <- kmeans(V_type2, centers = K, iter.max = 100, nstart = 20)
  
  # within cluster sum of sq
  wss_type1 <- k_means_output_type1$tot.withinss
  wss_type2 <- k_means_output_type2$tot.withinss
  return(list(k_means_output_type1,k_means_output_type2,wss_type1,wss_type2))
}

## not needed ###
# Function to compute misclustering error rate
misclustering_error_rate <- function(true_labels, predicted_labels) {
  # Construct the confusion matrix
  confusion_matrix <- table(true_labels, predicted_labels)
  
  # Solve the linear sum assignment problem (LSAP) to find the best permutation
  # of predicted clusters that minimizes the error.
  best_perm <- solve_LSAP(confusion_matrix, maximum = TRUE)
  
  # Rearrange the predicted labels based on the optimal permutation
  permuted_labels <- best_perm[predicted_labels]
  
  # Compute the number of misclassified points
  misclassified <- sum(true_labels != permuted_labels)
  
  # Calculate the misclustering error rate
  misclustering_rate <- misclassified / length(true_labels)
  
  return(misclustering_rate)
}


pMatrix.min <- function(A, B) { 
  n <- nrow(A) 
  D <- matrix(NA, n, n) 
  for (i in 1:n) { 
    for (j in 1:n) { 
      D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
    } } 
  vec <- c(solve_LSAP(D)) 
  list(A=A[vec,], pvec=vec) 
}

## hardcoded for 2 types
## i.e. T = 2
misclustering_error_rate_het<-function(node_info, A, T = 2, K, n){
  
  node_info_type1 = data.frame(node_info[1:(n/2),])
  colnames(node_info_type1) = c("type","block")
  true_labels_type1 = c(node_info_type1$block)
  
  node_info_type2 = data.frame(node_info[((n/2)+1):n,])
  colnames(node_info_type2) = c("type","block")
  true_labels_type2 = c(node_info_type2$block)
  
  ## Compute the regularized graph Laplacian L_tau
  L_tau <- gen_regularized_laplacian(A)
  
  V <- RSpectra::eigs_sym(L_tau, T*K, "LM")$vectors
  
  V_type1 = V[1:(n/2),]
  V_type2 = V[((n/2)+1):n,]
  
  # Apply k-means clustering to the eigenvectors
  k_means_output_type1 <- kmeans(V_type1, centers = K, iter.max = 1000, nstart = 100)
  k_means_output_type2 <- kmeans(V_type2, centers = K, iter.max = 1000, nstart = 100)
  
  predicted_label_type1 = k_means_output_type1$cluster
  predicted_label_type2 = k_means_output_type2$cluster
  
  # misclassified_type1 = misclustering_error_rate(true_labels_type1,predicted_label_type1)
  # misclassified_type2 = misclustering_error_rate(true_labels_type2,predicted_label_type2)
  
  ### type 1 ####
  E1 <- factor(predicted_label_type1, levels = 1:K)
  A <- table(E1, true_labels_type1)
  B <- diag(as.numeric(table(true_labels_type1)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type1 = sum(A) - sum(diag(A))
  
  ### type 2 ####
  E2 <- factor(predicted_label_type2, levels = 1:K)
  A <- table(E2, true_labels_type2)
  B <- diag(as.numeric(table(true_labels_type2)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type2 = sum(A) - sum(diag(A))
  
  
  return(list(number_of_misclassified_type1,number_of_misclassified_type2))
}


find_optimal_alpha <-function(A,X,true_node_labels,K){
  L_tau = gen_regularized_laplacian(A)
  
  # Compute eigenvalues
  eig_L <- eigen(tcrossprod(L_tau), symmetric = TRUE, only.values = TRUE)$values
  eig_X <- eigen(tcrossprod(X), symmetric = TRUE, only.values = TRUE)$values
  
  # Sort eigenvalues in decreasing order
  eig_L <- sort(eig_L, decreasing = TRUE)
  eig_X <- sort(eig_X, decreasing = TRUE)
  
  # Compute alpha_min
  alpha_min <- (eig_L[K] - eig_L[K + 1]) / eig_X[1]
  
  R = ncol(X)
  
  # Compute denominator depending on R and K
  if (R <= K) {
    denom <- eig_X[R]
  } else {
    denom <- eig_X[K] - eig_X[K + 1]
  }
  
  # Compute alpha_max
  alpha_max <- eig_L[1] / denom
  
  
  alpha_list <- seq(alpha_min, alpha_max, length.out = 10)  # 10 grid points
  no_of_misclassified_list <- numeric(length(alpha_list))
  
  for (i in seq_along(alpha_list)) {
    alpha <- alpha_list[i]
    L_tilde_temp = crossprod(L_tau) + alpha * tcrossprod(X)
    V = RSpectra::eigs_sym(L_tilde_temp, K, "LM")$vectors
    pred_node_labels = kmeans(V, centers = K, iter.max = 100, nstart = 20)$cluster
    
    E <- factor(pred_node_labels, levels = 1:K)
    A <- table(E, true_node_labels)
    B <- diag(as.numeric(table(true_node_labels)))
    A <- pMatrix.min(A, B)$A
    number_of_misclassified_temp = sum(A) - sum(diag(A))
    # error_temp = misclustering_error_rate(true_node_labels, pred_node_labels)
    no_of_misclassified_list[i] <- number_of_misclassified_temp
  }
  
  # Find alpha with minimum error
  alpha_optimal <- alpha_list[which.min(no_of_misclassified_list)]
  return(alpha_optimal)
}

## no of types = 2 is hardcoded here i.e we derive error rate for 2 types of nodes!
## this function gives errors for two types 
misclustering_error_changing_with_alpha <-function(alpha,node_info,L_tau_sq,X,T,K,n){
  
  node_info_type1 = data.frame(node_info[1:(n/2),])
  colnames(node_info_type1) = c("type","block")
  true_labels_type1 = c(node_info_type1$block)
  
  node_info_type2 = data.frame(node_info[(n/2+1):n,])
  colnames(node_info_type2) = c("type","block")
  true_labels_type2 = c(node_info_type2$block)
  
  ##### perform spectral clustering #####
  L_tilde = L_tau_sq + alpha * tcrossprod(X)
  V <- RSpectra::eigs_sym(L_tilde, T*K, "LM")$vectors
  
  V_type1 = V[1:(n/2),]
  V_type2 = V[(n/2+1):n,]
  
  # Apply k-means clustering to the eigenvectors
  k_means_output_type1 <- kmeans(V_type1, centers = K, iter.max = 1000, nstart = 100)
  k_means_output_type2 <- kmeans(V_type2, centers = K, iter.max = 1000, nstart = 100)
  
  predicted_label_type1 = k_means_output_type1$cluster
  predicted_label_type2 = k_means_output_type2$cluster
  
  # error_1 = misclustering_error_rate(true_labels_type1,predicted_label_type1)
  # error_2 = misclustering_error_rate(true_labels_type2,predicted_label_type2)
  
  ### type 1 ####
  E1 <- factor(predicted_label_type1, levels = 1:K)
  A <- table(E1, true_labels_type1)
  B <- diag(as.numeric(table(true_labels_type1)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type1 = sum(A) - sum(diag(A))
  
  ### type 2 ####
  E2 <- factor(predicted_label_type2, levels = 1:K)
  A <- table(E2, true_labels_type2)
  B <- diag(as.numeric(table(true_labels_type2)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type2 = sum(A) - sum(diag(A))
  
  
  return(list(number_of_misclassified_type1,number_of_misclassified_type2))
  
}

homogeneous_sbm <- function(r, n, K = 3){
  
  ones_vector_K <- rep(1, K)  # This creates a vector (1,1,...,1) of length K
  one_matrix_K <- ones_vector_K %*% t(ones_vector_K) ## K X K matrix with all ones
  
  ### probability controlling parameters###
  p1 = 0.2
  p2 = 0.2
  r1 = r
  r2 = r # = 0.2
  
  # Define the block probability matrices P11, P12, P21 (P21 = P12), P22
  P11 <- as.matrix(p1*one_matrix_K + r1*diag(1,K), K,K) # KxK matrix for type-1 to type-1 links
  #P12 <- as.matrix(p3*one_matrix_K + r3*diag(1,K), K,K)  # KxK matrix for type-1 to type-2 links
  P22 <- as.matrix(p2*one_matrix_K + r2*diag(1,K) ,K,K)  # KxK matrix for type-2 to type-2 links
  
  #########################################
  #set.seed(1234)
  # Assume equal partitioning of nodes into 6 type-block combinations for each type
  Z_type1 <- matrix(0, n/2, K)
  Z_type2 <- matrix(0, n/2, K)
  
  # Randomly assign each node to one of the 3 blocks for type-1
  for (i in 1:nrow(Z_type1)) {
    Z_type1[i, sample(1:K, 1)] <- 1
  }
  
  # Randomly assign each node to one of the 3 blocks for type-2
  for (i in 1:nrow(Z_type2)) {
    Z_type2[i, sample(1:K, 1)] <- 1
  }
  
  A1 = gen_adj(Z_type1,P11)
  A2 = gen_adj(Z_type2,P22)
  ##############################################################
  ## true labels ###
  true_labels_type1 = c()
  true_labels_type2 = c()
  for (i in 1:nrow(Z_type1)){
    true_labels_type1[i] = which(Z_type1[i,] == 1)
  }
  for (i in 1:nrow(Z_type2)){
    true_labels_type2[i] = which(Z_type2[i,] == 1)
  }
  
  L_tau1 = gen_regularized_laplacian(A1)
  L_tau2 = gen_regularized_laplacian(A2)
  
  V_type1 = get_eigenvector_matrix(L_tau1,K)
  V_type2 = get_eigenvector_matrix(L_tau2,K)
  
  predicted_label_type1 =  kmeans(V_type1, centers = K, iter.max = 100, nstart = 20)$cluster
  predicted_label_type2 =  kmeans(V_type2, centers = K, iter.max = 100, nstart = 20)$cluster
  
  ### type 1 ####
  E1 <- factor(predicted_label_type1, levels = 1:K)
  A <- table(E1, true_labels_type1)
  B <- diag(as.numeric(table(true_labels_type1)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type1 = sum(A) - sum(diag(A))
  
  ### type 2 ####
  E2 <- factor(predicted_label_type2, levels = 1:K)
  A <- table(E2, true_labels_type2)
  B <- diag(as.numeric(table(true_labels_type2)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type2 = sum(A) - sum(diag(A))
  
  
  return(list(number_of_misclassified_type1,number_of_misclassified_type2))
  
  # return(list(err1,err2))
  
}


homogeneous_cov_sbm<-function(d, r, n, K=3, T=2, R=3){
  
  ones_vector_K <- rep(1, K)  # This creates a vector (1,1,...,1) of length K
  one_matrix_K <- ones_vector_K %*% t(ones_vector_K) ## K X K matrix with all ones
  
  ### probability controlling parameters###
  p1 = 0.2
  p2 = 0.2
  r1 = r 
  r2 = r # = 0.2
  
  #r3 = 0.08
  #p3 = 0.5
  
  # Define the block probability matrices P11, P12, P21 (P21 = P12), P22
  P11 <- as.matrix(p1*one_matrix_K + r1*diag(1,K), K,K) # K x K matrix for type-1 to type-1 links
  #P12 <- as.matrix(p3*one_matrix_K + r3*diag(1,K), K,K)  # K x K matrix for type-1 to type-2 links
  P22 <- as.matrix(p2*one_matrix_K + r2*diag(1,K) ,K,K)  # K x K matrix for type-2 to type-2 links
  
  #########################################
  #set.seed(1234)
  # Assume equal partitioning of nodes into 6 type-block combinations for each type
  Z_type1 <- matrix(0, n/2, K)
  Z_type2 <- matrix(0, n/2, K)
  
  # Randomly assign each node to one of the 3 blocks for type-1
  for (i in 1:nrow(Z_type1)) {
    Z_type1[i, sample(1:K, 1)] <- 1
  }
  
  # Randomly assign each node to one of the 3 blocks for type-2
  for (i in 1:nrow(Z_type2)) {
    Z_type2[i, sample(1:K, 1)] <- 1
  }
  
  A1 = gen_adj(Z_type1,P11)
  A2 = gen_adj(Z_type2,P22)
  ##############################################################
  ## true labels ###
  true_labels_type1 = c()
  true_labels_type2 = c()
  for (i in 1:nrow(Z_type1)){
    true_labels_type1[i] = which(Z_type1[i,] == 1)
  }
  for (i in 1:nrow(Z_type2)){
    true_labels_type2[i] = which(Z_type2[i,] == 1)
  }
  
  L_tau1 = gen_regularized_laplacian(A1)
  L_tau2 = gen_regularized_laplacian(A2)
  
  ### covariate ######
  R = 3 ## dimension of each covariate 
  
  ##CAUION! dont set p_out_cov_type1 and p_out_cov_type2 greater than 0.1
  ##because will varry d from 0 upto 0.9.
  ## in that case p_in_cov would become bigger than 1!
  
  #### means for different combinations##
  p_out_cov_type1 = 0.08
  #p_in_cov_type1 = 0.2 + d
  
  p_out_cov_type2 = 0.05
  #p_in_cov_type1 = 0.3 + d
  
  M_type1 <- as.matrix(p_out_cov_type1*one_matrix_K + d*diag(1,K),nrow=T*K/2, ncol=R)
  M_type2 <- as.matrix(p_out_cov_type2*one_matrix_K + d*diag(1,K),nrow=T*K/2, ncol=R)
  
  ### expected value of X given Z (E(X|Z))
  expected_X_type1 = Z_type1 %*% M_type1
  expected_X_type2 = Z_type2 %*% M_type2
  
  ### generating the covariate vector at each node and covariate matrix is named as X
  X1 = matrix(rbinom(nrow(Z_type1) * R, 1, expected_X_type1), nrow = nrow(Z_type1), ncol = R)
  X2 = matrix(rbinom(nrow(Z_type2) * R, 1, expected_X_type2), nrow = nrow(Z_type2), ncol = R)
  ## no need to scale bernoulli variables
  
  optimal_alpha_type1 = find_optimal_alpha(A1,X1,true_labels_type1,K)
  optimal_alpha_type2 = find_optimal_alpha(A2,X2,true_labels_type2,K)
  
  ## prepare cov-assisted laplacians###
  L_tilde1 = crossprod(L_tau1) + optimal_alpha_type1*tcrossprod(X1)
  L_tilde2 = crossprod(L_tau2) + optimal_alpha_type2*tcrossprod(X2)
  
  V1 = RSpectra::eigs_sym(L_tilde1, K, which = "LM")$vectors
  V2 = RSpectra::eigs_sym(L_tilde2, K, which = "LM")$vectors
  
  pred_labels_type1 = kmeans(V1, centers = K, iter.max = 100, nstart = 20)$cluster
  pred_labels_type2 = kmeans(V2, centers = K, iter.max = 100, nstart = 20)$cluster
  
  
  ### type 1 ####
  E1 <- factor(pred_labels_type1, levels = 1:K)
  A <- table(E1, true_labels_type1)
  B <- diag(as.numeric(table(true_labels_type1)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type1 = sum(A) - sum(diag(A))
  
  ### type 2 ####
  E2 <- factor(pred_labels_type2, levels = 1:K)
  A <- table(E2, true_labels_type2)
  B <- diag(as.numeric(table(true_labels_type2)))
  A <- pMatrix.min(A, B)$A
  number_of_misclassified_type2 = sum(A) - sum(diag(A))
  
  
  return(list(number_of_misclassified_type1,number_of_misclassified_type2))
  
  # err1 = misclustering_error_rate(true_labels_type1,pred_labels_type1)
  # err2 = misclustering_error_rate(true_labels_type2,pred_labels_type2)
  # return(list(err1,err2))
  
}

heterogeneous_sbm <- function(r, r3, n, K=3, T=2){
  #n = 250 # no of nodes
  #K = 3 # no of blocks
  #T = 2 # no of types of blocks, i.e. two types of blocks- type1,type2
  
  ones_vector_K <- rep(1, K)  # This creates a vector (1,1,...,1) of length K
  
  one_matrix_K <- ones_vector_K %*% t(ones_vector_K) ## K X K matrix with all ones
  
  ### probability controlling parameters###
  p1 = 0.2
  p2 = 0.2
  r1 = r
  r2 = r
  # r3 = 0.08
  #### vary r3 accordingly upto 0.85 so that p3+r3 may not exist 1
  p3 = 0.15
  
  # Define the block probability matrices P11, P12, P21 (P21 = P12), P22
  P11 <- as.matrix(p1*one_matrix_K + r1*diag(1,K), K,K) # KxK matrix for type-1 to type-1 links
  P12 <- as.matrix(p3*one_matrix_K + r3*diag(1,K), K,K)  # KxK matrix for type-1 to type-2 links
  P22 <- as.matrix(p2*one_matrix_K + r2*diag(1,K) ,K,K)  # KxK matrix for type-2 to type-2 links
  
  # Construct the full probability matrix P
  
  P <- rbind(cbind(P11, P12),
             cbind(t(P12), P22))
  
  #########################################
  #set.seed(1234)
  # Assume equal partitioning of nodes into 6 type-block combinations for each type
  Z_type1 <- matrix(0, n/2, K)
  Z_type2 <- matrix(0, n/2, K)
  
  # Randomly assign each node to one of the 3 blocks for type-1
  for (i in 1:nrow(Z_type1)) {
    Z_type1[i, sample(1:(K), 1)] <- 1
  }
  
  # Randomly assign each node to one of the 3 blocks for type-2
  for (i in 1:nrow(Z_type2)) {
    Z_type2[i, sample(1:(K), 1)] <- 1
  }
  zero_matrix = matrix(0,n/2,K)
  
  Z_type1 = cbind(Z_type1,zero_matrix)
  Z_type2 = cbind(zero_matrix, Z_type2)
  # Combine type-1 and type-2 membership matrices
  Z <- rbind(Z_type1, Z_type2)
  #colSums(Z) ### to check how many nodes in each (type-block) combinations
  ##############################################
  # Generate the adjacency matrix A from expected_A (E(A|Z)) as parameter
  A = gen_adj(Z,P)
  
  ####################################################################
  # Create a vector to store ground truth node information(node type and block membership)
  node_info <- matrix(0, n, 2)  # Column 1: type, Column 2: block
  
  # Assign type and block information
  for (i in 1:n) {
    ## type membership
    node_info[i, 1] <- ifelse(which(Z[i,]==1)<=(T*K/2), 1, 2)  # first 3 positions of a row of Z:type1; last 3 positions of that row:type2
    ## Block membership
    if (which(Z[i,] == 1) %in% c(1, 4)) {
      
      node_info[i,2] <- 1  # Block 1
      
    } else if (which(Z[i,] == 1) %in% c(2, 5)) {
      
      node_info[i,2] <- 2  #  Block 2
      
    } else if (which(Z[i,] == 1) %in% c(3, 6)) {
      
      node_info[i,2] <- 3  # Block 3
      
    }
  }
  
  ### yes_het_no_cov_error #####
  
  # error_type1_het = misclustering_error_rate_het(node_info = node_info, A = A, T = 2, K=K, n = n)[[1]]
  # error_type2_het = misclustering_error_rate_het(node_info = node_info, A = A, T = 2, K=K, n = n)[[2]]
  # 
  # return(list(error_type1_het,error_type2_het))
  number_of_misclassified_type1_het = misclustering_error_rate_het(node_info = node_info, A = A, T = 2, K=K, n = n)[[1]]
  number_of_misclassified_type2_het = misclustering_error_rate_het(node_info = node_info, A = A, T = 2, K=K, n = n)[[2]]
  
  return(list(number_of_misclassified_type1_het,number_of_misclassified_type2_het))
  
}

het_cov_sbm <-function(d, r, r3, n, K=3, T=2, R=3){
  #n = 250 # no of nodes
  #K = 3 # no of blocks
  #T = 2 # no of types of blocks, i.e. two types of blocks- type1,type2
  
  ones_vector_K <- rep(1, K)  # This creates a vector (1,1,...,1) of length K
  
  one_matrix_K <- ones_vector_K %*% t(ones_vector_K) ## K X K matrix with all ones
  
  ### probability controlling parameters###
  p1 = 0.2
  p2 = 0.2
  r1 = r
  r2 = r
  # r3 = 0.08
  #### vary r3 accordingly upto 0.85 so that p3+r3 may not exist 1
  p3 = 0.15
  
  # Define the block probability matrices P11, P12, P21 (P21 = P12), P22
  P11 <- as.matrix(p1*one_matrix_K + r1*diag(1,K), K,K) # KxK matrix for type-1 to type-1 links
  P12 <- as.matrix(p3*one_matrix_K + r3*diag(1,K), K,K)  # KxK matrix for type-1 to type-2 links
  P22 <- as.matrix(p2*one_matrix_K + r2*diag(1,K) ,K,K)  # KxK matrix for type-2 to type-2 links
  
  # Construct the full probability matrix P
  
  P <- rbind(cbind(P11, P12),
             cbind(t(P12), P22))
  
  #########################################
  #set.seed(1234)
  # Assume equal partitioning of nodes into 6 type-block combinations for each type
  Z_type1 <- matrix(0, n/2, K)
  Z_type2 <- matrix(0, n/2, K)
  
  # Randomly assign each node to one of the 3 blocks for type-1
  for (i in 1:nrow(Z_type1)) {
    Z_type1[i, sample(1:(K), 1)] <- 1
  }
  
  # Randomly assign each node to one of the 3 blocks for type-2
  for (i in 1:nrow(Z_type2)) {
    Z_type2[i, sample(1:(K), 1)] <- 1
  }
  zero_matrix = matrix(0,n/2,K)
  
  Z_type1 = cbind(Z_type1,zero_matrix)
  Z_type2 = cbind(zero_matrix, Z_type2)
  # Combine type-1 and type-2 membership matrices
  Z <- rbind(Z_type1, Z_type2)
  #colSums(Z) ### to check how many nodes in each (type-block) combinations
  ##############################################
  
  # Generate the adjacency matrix A from expected_A (E(A|Z)) as parameter
  A = gen_adj(Z,P)
  
  ####################################################################
  # Create a vector to store ground truth node information(node type and block membership)
  node_info <- matrix(0, n, 2)  # Column 1: type, Column 2: block
  
  # Assign type and block information
  for (i in 1:n) {
    ## type membership
    node_info[i, 1] <- ifelse(which(Z[i,]==1)<=(T*K/2), 1, 2)  # first 3 positions of a row of Z:type1; last 3 positions of that row:type2
    ## Block membership
    if (which(Z[i,] == 1) %in% c(1, 4)) {
      
      node_info[i,2] <- 1  # Block 1
      
    } else if (which(Z[i,] == 1) %in% c(2, 5)) {
      
      node_info[i,2] <- 2  #  Block 2
      
    } else if (which(Z[i,] == 1) %in% c(3, 6)) {
      
      node_info[i,2] <- 3  # Block 3
      
    }
  }
  
  #############################################################################
  ### covariate ######
  R = 3 ## dimension of each covariate 
  
  ##CAUION! dont set p_out_cov_type1 and p_out_cov_type2 greater than 0.1
  ##because you will varry d from 0 upto 0.9.
  ## in that case p_in_cov would become bigger than 1!
  #### means for different combinations##
  p_out_cov_type1 = 0.08
  #p_in_cov_type1 = 0.08 + d
  
  p_out_cov_type2 = 0.05
  #p_in_cov_type1 = 0.05 + d
  
  M_type1 <- as.matrix(p_out_cov_type1*one_matrix_K + d*diag(1,K),nrow=T*K/2, ncol=R)
  M_type2 <- as.matrix(p_out_cov_type2*one_matrix_K + d*diag(1,K),nrow=T*K/2, ncol=R)
  
  M = rbind(M_type1,M_type2)
  
  #print(M)
  ### expected value of X given Z (E(X|Z))
  expected_X = Z%*%M
  
  ### generating the covariate vector at each node and covariate matrix is named as X
  X = matrix(rbinom(n * R, 1, expected_X), nrow = n, ncol = R)
  ## no need to scale bernoulli variables
  #X = scale(X)
  
  ### calculating degrees of the nodes
  degrees = rowSums(A)
  
  ## calculating the regularizer tau as average node degree
  tau = mean(degrees)
  
  ### regularized degree matrix
  #D_tau = diag(degrees+tau)  
  
  ## D_tau^{-1/2}
  D_tau_power_neg_half = diag(1/sqrt(degrees+tau))
  
  ## Compute the regularized graph Laplacian L_tau
  L_tau <-  D_tau_power_neg_half %*% A %*% D_tau_power_neg_half
  #################################################
  #### calculate L_tau_sq and X_sq
  L_tau_sq = crossprod(L_tau)
  XX_t = tcrossprod(X)
  
  ###################################################
  ####calculate tuning parameter alpha########## 
  #Compute the eigenvalues 
  eigenvalues_of_L_tau_sq <- eigen(L_tau_sq)$values
  eigenvalues_of_XX_t <- eigen(XX_t)$values
  ##Sort the eigenvalues in decreasing order
  sorted_eigenvalues_of_L_tau_sq <- sort(eigenvalues_of_L_tau_sq, decreasing = TRUE)
  sorted_eigenvalues_of_XX_t <- sort(eigenvalues_of_XX_t, decreasing = TRUE)
  #################################
  num_alpha_min = (sorted_eigenvalues_of_L_tau_sq[T*K] - sorted_eigenvalues_of_L_tau_sq[T*K+1])
  denom_alpha_min = sorted_eigenvalues_of_XX_t[1]
  alpha_min = num_alpha_min/denom_alpha_min
  
  num_alpha_max = sorted_eigenvalues_of_L_tau_sq[1]
  denom_alpha_max = sorted_eigenvalues_of_XX_t[R] * ifelse(R<=T*K,1,0) + (sorted_eigenvalues_of_XX_t[T*K]-sorted_eigenvalues_of_XX_t[T*K+1])*ifelse(R>T*K,1,0)
  alpha_max = num_alpha_max/denom_alpha_max
  #####################################################################
  
  ### finding optimal alpha
  range = alpha_max - alpha_min
  # xval = seq(alpha_min, alpha_max, range/50)
  # yval_1 = sapply(xval, function(x) wss_alpha(alpha=x, L_tau_sq = L_tau_sq , X = X, T =T, K=K, n=n)[[3]])
  # yval_2 = sapply(xval, function(x) wss_alpha(alpha=x, L_tau_sq = L_tau_sq , X = X, T =T, K=K, n=n)[[4]])
  
  ############### misclustering error############
  alpha_vals =  seq(alpha_min, alpha_max, range/5)
  
  ## yes_het_yes_cov error ##
  number_of_misclassified_type1 = sapply(alpha_vals,function(x) misclustering_error_changing_with_alpha(alpha = x, node_info = node_info, L_tau_sq =  L_tau_sq, X = X, T = T, K=K, n=n)[[1]])
  number_of_misclassified_type2 = sapply(alpha_vals,function(x) misclustering_error_changing_with_alpha(alpha = x, node_info = node_info, L_tau_sq =  L_tau_sq, X = X, T = T, K=K, n=n)[[2]])
  
  
  return(list(min(number_of_misclassified_type1),min(number_of_misclassified_type2)))
  
}

