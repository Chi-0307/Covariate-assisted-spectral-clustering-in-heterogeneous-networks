###############################################################################################################
###############################################################################################################
#=== Descriptions =======
## contains data prep for unga data 1990-2000 
## along with hom,het,homcov, hetcov clustering of the data
## ========================================================
library(unvotes)
options(timeout = 300)   # seconds (5 minutes)
library(WDI)
library(data.table);library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Matrix)
library(ggrepel)
library(Matrix)
library(RSpectra)
library(mclust)## ARI Calculation
library(maps)



####### HET_COV + HET_NO_COV ##########################
get_wdi_colnames <- function(indicators) {
  
  # load these 
  if (!requireNamespace("WDI", quietly = TRUE))
    stop("Package 'WDI' is required")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    stop("Package 'stringr' is required")
  
  # get metadata
  WDIcache <- WDI::WDIcache()
  ind_meta <- WDIcache$series
  
  # helper to clean names
  make_colname <- function(x) {
    x |>
      tolower() |>
      stringr::str_replace_all("[^a-z0-9]+", "_") |>
      stringr::str_replace_all("^_|_$", "")
  }
  
  # map indicator -> clean column name
  out <- ind_meta[ind_meta$indicator %in% indicators, c("indicator", "name")]
  out$colname <- make_colname(out$name)
  
  # keep order same as input
  out <- out[match(indicators, out$indicator), ]
  
  out$colname
}

## ================= ONE-TIME SETUP (run once) =================================
## because when the function is run in parallel it hits WDI-API multiple times, so timeout happens
##  indicators given 
## ===========================================================
years <- 1990:2000

wdi_ind_list <- c(
  "NY.GDP.PCAP.CD",
  "NE.TRD.GNFS.ZS",
  "BX.KLT.DINV.WD.GD.ZS"
)

wdi_all <- WDI(
  country   = "all",
  indicator = wdi_ind_list,
  start     = min(years),
  end       = max(years),
  extra     = TRUE,
  cache     = NULL
)
saveRDS(wdi_all, "data/wdi_all_years_all_indicators.rds")

################################################

wdi_colnames <- get_wdi_colnames(wdi_ind_list)

saveRDS(wdi_colnames,
        "data/wdi_colnames.rds")

## ============= PRE-LOAD (ONCE PER SESSION) =====================

##### paths #######
### save it to use inside the function ####
wdi_rds_path  <- "data/wdi_all_years_all_indicators.rds"
wdi_col_path  <- "data/wdi_colnames.rds"
hdi_csv_path  <- "data/HDR25_Composite_indices_complete_time_series.csv"

# load data
wdi_data_all       <- readRDS(wdi_rds_path)
wdi_colnames_map   <- readRDS(wdi_col_path)
hdi_data_complete  <- fread(hdi_csv_path)

# name the column map by indicator
wdi_ind_list <- c(
  "NY.GDP.PCAP.CD",
  "NE.TRD.GNFS.ZS",
  "BX.KLT.DINV.WD.GD.ZS"
)
names(wdi_colnames_map) <- wdi_ind_list



#===============================================================================
##======== function for het-cov and het-only clustering ========================
## =============================================================================
#### method = het_yes_cov gives het_cov function; method = het_no_cov gives het function 

het_cov_unga_final <-function(input_year, K, method = c('het_yes_cov','het_no_cov'), wdi_indicators, include_hdi){
  
  ############################################################################################
  ####################### prepare data #######################################################
  ############################################################################################
  
  ###### wdi data and hdi data prep ############
  
  wdi_data <- wdi_data_all
  wdi_data <- wdi_data[, !names(wdi_data) %in% "lastupdated"]

  
  ## in un_votes the name is Laos,Brunei,Myanmar (Burma), Vietnam resp. Need to match for joining. So we changed the name here.
  wdi_data$country[wdi_data$country == "Lao PDR"] <- "Laos"
  wdi_data$country[wdi_data$country == "Brunei Darussalam"] <- "Brunei"
  wdi_data$country[wdi_data$country == "Myanmar"] <- "Myanmar (Burma)"
  wdi_data$country[wdi_data$country == "Viet Nam"] <- "Vietnam"
  
  
  wdi_data <- subset(wdi_data, select = -c(capital,longitude,latitude,income,lending,status))
  
  # subset year + required cols
  wdi_data <- subset(
    wdi_data,
    year == input_year,
    select = c("country","iso2c","iso3c","year","region", wdi_indicators)
  )
  
  # rename indicators : clean names (offline)
  wdi_cov <- unname(wdi_colnames_map[wdi_indicators])
  colnames(wdi_data)[colnames(wdi_data) %in% wdi_indicators] <- wdi_cov
  
  ### keep wdi_cov columns along with few other necessary cols 

  wdi_data_selected <- subset(wdi_data, !region %in% c("Aggregates")) ## removed because 'Aggregates' rows are regions not countries for ex sub-saharan africa
  
  if (include_hdi == 'yes'){
    
    # to join with wdi both df shud have the same iso3 name
    names(hdi_data_complete)[names(hdi_data_complete) == "iso3"] <- "iso3c"
    
    hdi_col <- paste0("hdi_", input_year);gdi_col<-paste0("gdi_", input_year);
    
    # subset necessary columns for input year

    necessary_cols = c("iso3c","country","region", hdi_col, gdi_col)
    hdi_data_selected <- subset(hdi_data_complete, select = necessary_cols)
    
    wdi_hdi_joined <- inner_join(wdi_data_selected, hdi_data_selected, by = "iso3c")
    to_join_with_un <- wdi_hdi_joined
    
  }else{
    
    to_join_with_un <- wdi_data_selected
    
  }
  
  #########################################################################
  #########################################################################
  
  ##### un data prep #############
  
  un_roll_calls = data.frame(un_roll_calls)
  un_roll_call_issues = data.frame(un_roll_call_issues)
  un_votes = data.frame(un_votes)
  
  d1 = full_join(un_roll_calls, un_roll_call_issues, by = "rcid")
  d2 = inner_join(d1, un_votes, by = "rcid",relationship = "many-to-many")

  
  d2_selected = d2[year(d2$date) == input_year,]
  names(d2_selected)[names(d2_selected) == "country_code"] <- "iso2c"

  
  ##### join wdi+hdi(if taken)+un_vote
  ddd = inner_join(d2_selected, to_join_with_un, by = 'iso2c')
  
  ddd_issues_imputed <- ddd %>%
    mutate(
      issue = case_when(
        !is.na(issue) ~ issue,  # keep existing values
        
        grepl("Palestin|Syrian Golan|June 1967|occupied", descr, ignore.case = TRUE) ~ "Palestinian conflict",
        grepl("human cloning|discrimination|women|democratization|defamation of religions|displaced", descr, ignore.case = TRUE) ~ "Human rights",
        grepl("economic|development|deserts|desertification", descr, ignore.case = TRUE) ~ "Economic development",
        grepl("colonialism|decolonization", descr, ignore.case = TRUE) ~ "Colonialism",
        grepl("armaments|outer space|confidence-building", descr, ignore.case = TRUE) ~ "Arms control and disarmament",
        grepl("atomic energy|nuclear", descr, ignore.case = TRUE) ~ "Nuclear weapons and nuclear material",
        grepl("oceans|environment", descr, ignore.case = TRUE) ~ "Environment",
        
        is.na(issue) ~ "Other" # fallback if nothing matches
      )
    )
  
  
  if (include_hdi == 'yes'){
    my_df <- ddd_issues_imputed %>%
      
      rename(
        hdi = paste0("hdi_", input_year),
        gdi = paste0("gdi_", input_year)
      )
    
    df_imputed <- my_df %>%
      group_by(region.x) %>%
      mutate(
        across(
          all_of(wdi_cov),
          ~ ifelse(
            is.na(.),
            round(mean(., na.rm = TRUE), 3),
            .
          )
        ),
        hdi = ifelse(is.na(hdi), round(mean(hdi, na.rm = TRUE),3), hdi),
        gdi = ifelse(is.na(gdi), round(mean(hdi, na.rm = TRUE),3), gdi)
      ) %>%
      ungroup()
  }else{
    
    my_df = ddd_issues_imputed
    
    df_imputed <- my_df %>%
      group_by(region) %>%
      mutate(
        across(
          all_of(wdi_cov),
          ~ ifelse(
            is.na(.),
            round(mean(., na.rm = TRUE), 3),
            .
          )
        )
      ) %>%
      ungroup()
    
  }
  
  
  ##################################################################
  ############### prepare X ########################################
  ##################################################################
  
  #Get list of countries and resolutions
  if(include_hdi == 'yes'){
    countries <- sort(unique(df_imputed$country))
    
  }else{
    countries <- sort(unique(df_imputed$country.x))
    names(df_imputed)[names(df_imputed) == 'country.x'] <- 'country'
  }
  
  descriptions <- sort(unique(df_imputed$descr))
  
  ### cleaning descriptions and df_imputed$descr as they may contain " " or NA...
  descriptions <- descriptions[descriptions != "" & !is.na(descriptions)]
  df_imputed <- df_imputed[df_imputed$descr != "" & !is.na(df_imputed$descr), ]
  
  n_country <- length(countries)
  n_descr <- length(descriptions)
  
  if(include_hdi == 'yes'){
    cont_covariates = c(wdi_cov, "hdi","gdi")
  }else{
    cont_covariates = c(wdi_cov)
  }
  
  categorical_covariates = c(unique(df_imputed$issue))
  
  #Create country covariate matrix (n_country x 8)
  X_country <- matrix(0, nrow = n_country, ncol = length(cont_covariates) + length(categorical_covariates))
  
  all_cov = c(cont_covariates,categorical_covariates)
  
  for(i in 1:n_country){
    c <- countries[i]
    
    # Subset rows for this country
    sub_df <- df_imputed[df_imputed$country == c, ] 
    
    ### first fill wdi covariates
    X_country[i, 1:length(wdi_cov)] <-
      sapply(wdi_cov, function(v) {
        x <- sub_df[[v]]
        ifelse(any(!is.na(x)), x[!is.na(x)][1], 0)
      })
    
    ## then fill hdi and gdi if include_hdi is "yes"
    if (include_hdi == 'yes'){
      X_country[i, (length(wdi_cov) + 1)] <- sub_df$hdi[!is.na(sub_df$hdi)][1]
      X_country[i, (length(wdi_cov) + 2)] <- sub_df$gdi[!is.na(sub_df$gdi)][1]
    }
  }   
  
  #The last columns are 0 for countries
  
  # Define issue categories (5 types)
  issue_list <- sort(unique(na.omit(df_imputed$issue)))  # e.g. Environmental, Human Rights, etc.
  
  # Create resolution covariate matrix (n_descr x (#cont_cov+#cat_cov))
  X_descr <- matrix(0, nrow = n_descr, ncol = length(cont_covariates)+length(categorical_covariates))
  
  for(j in 1:n_descr){
    d <- descriptions[j]
    
    # Get all issues linked to this resolution
    issues_for_d <- unique(df_imputed$issue[df_imputed$descr == d])
    
    # Set 1 in corresponding positions
    for(iss in issues_for_d){
      if(!is.na(iss)){
        index <- match(iss, issue_list)  
        if(!is.na(index)){
          X_descr[j, length(cont_covariates) + index] <- 1
        }
      }
    }
  }
  
  # Stack the matrices
  X <- rbind(X_country, X_descr)
  
  
  # Assign column names
  if(include_hdi == 'yes'){
    colnames(X) <- c(wdi_cov, "hdi", "gdi", issue_list)
  }else{
    colnames(X) <- c(wdi_cov, issue_list)
  }
  
  # Assign row names to X
  # Clean resolution names if needed (optional: shorten)
  short_descr <- substr(descriptions, 1, 35)  
  
  # Combine country names and descr names as rownames
  rownames(X) <- c(countries, short_descr)
  
  # Create type vector
  type_vec <- c(rep("country", n_country), rep("descr", n_descr))
  
  scaling <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  # scale only continuous covariates
  X[1:n_country, cont_covariates] <- apply(
    X[1:n_country, cont_covariates, drop = FALSE],
    2,
    scaling
  )
  
  X[is.na(X)] <- 0
  
  
  #######################################################################
  ####################### prepare adj ###################################
  #######################################################################
  
  # Index maps
  country_index <- setNames(seq_along(countries), countries)
  descr_index <- setNames(seq_along(descriptions), descriptions)
  
  # Initialize binary bipartite matrix
  B <- matrix(0, nrow = n_country, ncol = n_descr)
  
  for (i in 1:nrow(df_imputed)) {
    vote <- df_imputed$vote[i]
    if (!is.na(vote) && vote == "yes") {
      c <- df_imputed$country[i]
      d <- df_imputed$descr[i]
      row <- country_index[[c]]
      col <- descr_index[[d]]
      B[row, col] <- 1
    }
  }
  # Build full symmetric adjacency matrix : bipartite(heterogeneous)graph with nodes: country and issues
  A <- matrix(0, nrow = (n_country + n_descr), ncol = (n_country + n_descr))
  A[1:n_country, (n_country + 1):(n_country + n_descr)] <- B
  A[(n_country + 1):(n_country + n_descr), 1:n_country] <- t(B)
  
  # assign row/col names
  rownames(A) <- c(countries, paste0("res_", seq_along(descriptions)))
  colnames(A) <- rownames(A)
  
  n = n_country + n_descr
  
  tau = mean(rowSums(A))
  
  inv_sqrt_deg <- ifelse(rowSums(A) > 0, 1 / sqrt(rowSums(A + tau * diag(n))), 0)
  
  D_inv_sqrt <- diag(inv_sqrt_deg)
  
  L_sym = tcrossprod(D_inv_sqrt %*% A, D_inv_sqrt) 
  R = ncol(X)
  
  #K = 3
  T = 2 ## no of types, country+resolution
  
  ### choose alpha this way ##
  LL    <- crossprod(L_sym)                       # L_tau L_tau  (= L^2, symmetric PSD)
  X_xt  <- tcrossprod(X); X_xt[is.na(X_xt)] <- 0   # XX^T
  
  lam1_LL <- eigs_sym(LL,   k = 1, which = "LM")$values[1]   # leading eigenvalue of L^2
  lam1_X  <- eigs_sym(X_xt, k = 1, which = "LM")$values[1]   # leading eigenvalue of XX^T
  
  alpha0  <- lam1_LL / lam1_X                      # Binkiewicz balancing value
  
  
  if(method== 'het_yes_cov'){
    # Compute XX^T
    X_xt <- tcrossprod(X)
    # Replace any NaN by 0
    X_xt[is.na(X_xt)] <- 0
    L_tilde = crossprod(L_sym) + alpha0 * X_xt 
    
    ## replace na by 0 
    L_tilde[is.na(L_tilde)] <- 0
    
    eig <- RSpectra::eigs_sym(L_tilde, k = T*K, which = "LM")  # "LM" = largest magnitude
    
    U <- eig$vectors  # n x TK matrix
    
    # Row-wise normalization
    row_norms <- sqrt(rowSums(U^2))
    row_norms[row_norms == 0] <- 1
    
    U_norm <- U / row_norms
    
    U_norm_country = U_norm[1:n_country,]

    # Run k-means on rows of U_country
    set.seed(3456) 
    
    # Output
    kmeans_result <- kmeans(U_norm_country, centers = K, nstart = 10)
    kmeans_output_cluster <- kmeans_result$cluster
    
  }else if(method == 'het_no_cov'){
    
    ### replace NA by 0
    L_sym[is.na(L_sym)] <- 0
    eig <- RSpectra::eigs_sym(L_sym, k = T*K, which = "LM")  # "LM" = largest magnitude
    
    U <- eig$vectors  # n x TK matrix
    
    #Row-wise normalize (optional) 
    # Row-wise normalization
    row_norms <- sqrt(rowSums(U^2))
    row_norms[row_norms == 0] <- 1
    U_norm <- U / row_norms
    
    U_norm_country = U_norm[1:n_country,]
    
    # Run k-means on rows of U_country
    set.seed(3456)
    
    # Output
    kmeans_result <- kmeans(U_norm_country, centers = K, nstart = 10)
    kmeans_output_cluster <- kmeans_result$cluster
  }
  
  country_clusters <- kmeans_output_cluster[1:n_country]
  
  # Map country names to clusters
  cluster_map <- data.frame(country = countries, cluster = country_clusters, 
                            country_degrees = colSums(A)[1:n_country])
  
  return(cluster_map)
  
}



#===============================================================================
##============= function for hom-cov and hom-only clustering ===================
## =============================================================================

############## HOM-COV ####################################################

### (similarity functions)##########
#cosine similarity #
normalize <- function(X){
  ## sometimes it may happen some countries give NV to all resolutions.
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

get_wdi_colnames <- function(indicators) {
  
  # load these 
  if (!requireNamespace("WDI", quietly = TRUE))
    stop("Package 'WDI' is required")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    stop("Package 'stringr' is required")
  
  # get metadata
  WDIcache <- WDI::WDIcache()
  ind_meta <- WDIcache$series
  
  # helper to clean names
  make_colname <- function(x) {
    x |>
      tolower() |>
      stringr::str_replace_all("[^a-z0-9]+", "_") |>
      stringr::str_replace_all("^_|_$", "")
  }
  
  # map indicator -> clean column name
  out <- ind_meta[ind_meta$indicator %in% indicators, c("indicator", "name")]
  out$colname <- make_colname(out$name)
  
  # keep order same as input
  out <- out[match(indicators, out$indicator), ]
  
  out$colname
}

## ================= ONE-TIME SETUP (run once) ------------------------
## because when the function is run in parallel it hits WDI-API multiple times, so timeout happens
##  indicators given 
##  =================
years <- 1990:2000

wdi_ind_list <- c(
  "NY.GDP.PCAP.CD",
  "NE.TRD.GNFS.ZS",
  "BX.KLT.DINV.WD.GD.ZS"
)

wdi_all <- WDI(
  country   = "all",
  indicator = wdi_ind_list,
  start     = min(years),
  end       = max(years),
  extra     = TRUE,
  cache     = NULL
)
saveRDS(wdi_all, "data/wdi_all_years_all_indicators.rds")

# ################################################

wdi_colnames <- get_wdi_colnames(wdi_ind_list)
saveRDS(wdi_colnames,
        "data/wdi_colnames.rds")

## ===================== PRE-LOAD (ONCE PER SESSION) =====================

# paths
### save it to use inside the function ####
wdi_rds_path  <- "data/wdi_all_years_all_indicators.rds"
wdi_col_path  <- "data/wdi_colnames.rds"
hdi_csv_path  <- "data/HDR25_Composite_indices_complete_time_series.csv"

# load data
wdi_data_all       <- readRDS(wdi_rds_path)
wdi_colnames_map   <- readRDS(wdi_col_path)
hdi_data_complete  <- fread(hdi_csv_path)

# name the column map by indicator
wdi_ind_list <- c(
  "NY.GDP.PCAP.CD",
  "NE.TRD.GNFS.ZS",
  "BX.KLT.DINV.WD.GD.ZS"
)
names(wdi_colnames_map) <- wdi_ind_list


###########################################################################
############ HOM-COV ######################################################
###########################################################################

hom_cov_unga_final <- function(input_year, similarity_measure, threshold = c('q1','q2','q3','mean'), K, wdi_indicators, include_hdi){
  
  ############## prepare data ###############
  ####################################################
  ## subsetting the wdi data for given input_year ####
  ####################################################
  
  ###### wdi data and hdi data prep ############
  
  wdi_data <- wdi_data_all
  wdi_data <- wdi_data[, !names(wdi_data) %in% "lastupdated"]
  
  ## in un_votes the name is Laos,Brunei,Myanmar (Burma),resp. Need to match for joining. So we changed the name here.
  wdi_data$country[wdi_data$country == "Lao PDR"] <- "Laos"
  wdi_data$country[wdi_data$country == "Brunei Darussalam"] <- "Brunei"
  wdi_data$country[wdi_data$country == "Myanmar"] <- "Myanmar (Burma)"
  wdi_data$country[wdi_data$country == "Viet Nam"] <- "Vietnam"
  
  
  wdi_data <- subset(wdi_data, select = -c(capital,longitude,latitude,income,lending,status))
  
  # subset year + required cols
  wdi_data <- subset(
    wdi_data,
    year == input_year,
    select = c("country","iso2c","iso3c","year","region", wdi_indicators)
  )
  
  # rename indicators : clean names (offline)
  wdi_cov <- unname(wdi_colnames_map[wdi_indicators])
  colnames(wdi_data)[colnames(wdi_data) %in% wdi_indicators] <- wdi_cov
  
  
  wdi_data_selected <- subset(wdi_data, !region %in% c("Aggregates")) ## removed because 'Aggregates' rows are regions not countries for ex sub-saharan africa
  
  if (include_hdi == 'yes'){
    
    # to join with wdi both df shud have the same iso3 name
    names(hdi_data_complete)[names(hdi_data_complete) == "iso3"] <- "iso3c"
    
    hdi_col <- paste0("hdi_", input_year);gdi_col<-paste0("gdi_", input_year);
      
    necessary_cols = c("iso3c","country","region", hdi_col, gdi_col)
    hdi_data_selected <- subset(hdi_data_complete, select = necessary_cols)
    
    wdi_hdi_joined <- inner_join(wdi_data_selected, hdi_data_selected, by = "iso3c")
    to_join_with_un <- wdi_hdi_joined
    
  }else{
    
    to_join_with_un <- wdi_data_selected
    
  }
  
  #########################################################################
  ############### un data prep ############################################
  
  un_roll_calls = data.frame(un_roll_calls)
  un_roll_call_issues = data.frame(un_roll_call_issues)
  un_votes = data.frame(un_votes)
  
  d1 = full_join(un_roll_calls, un_roll_call_issues, by = "rcid")
  d2 = inner_join(d1, un_votes, by = "rcid",relationship = "many-to-many")
  
  
  d2_selected = d2[year(d2$date) == input_year,]
  names(d2_selected)[names(d2_selected) == "country_code"] <- "iso2c"

  
  ##### join wdi+hdi(if taken)+un_vote
  ddd = inner_join(d2_selected, to_join_with_un, by = 'iso2c')
  
  ddd_issues_imputed <- ddd %>%
    mutate(
      issue = case_when(
        !is.na(issue) ~ issue,  # keep existing values
        
        grepl("Palestin|Syrian Golan|June 1967|occupied", descr, ignore.case = TRUE) ~ "Palestinian conflict",
        grepl("human cloning|discrimination|women|democratization|defamation of religions|displaced", descr, ignore.case = TRUE) ~ "Human rights",
        grepl("economic|development|deserts|desertification", descr, ignore.case = TRUE) ~ "Economic development",
        grepl("colonialism|decolonization", descr, ignore.case = TRUE) ~ "Colonialism",
        grepl("armaments|outer space|confidence-building", descr, ignore.case = TRUE) ~ "Arms control and disarmament",
        grepl("atomic energy|nuclear", descr, ignore.case = TRUE) ~ "Nuclear weapons and nuclear material",
        grepl("oceans|environment", descr, ignore.case = TRUE) ~ "Environment",
        
        is.na(issue) ~ "Other" # fallback if nothing matches
      )
    )
  
  
  if (include_hdi == 'yes'){
    my_df <- ddd_issues_imputed %>%
      
      rename(
        hdi = paste0("hdi_", input_year),
        gdi = paste0("gdi_", input_year)
      )
    
    df_imputed <- my_df %>%
      group_by(region.x) %>%
      mutate(
        across(
          all_of(wdi_cov),
          ~ ifelse(
            is.na(.),
            round(mean(., na.rm = TRUE), 3),
            .
          )
        ),
        hdi = ifelse(is.na(hdi), round(mean(hdi, na.rm = TRUE),3), hdi),
        gdi = ifelse(is.na(gdi), round(mean(hdi, na.rm = TRUE),3), gdi)
      ) %>%
      ungroup()
  }else{
    
    my_df = ddd_issues_imputed
    
    df_imputed <- my_df %>%
      group_by(region) %>%
      mutate(
        across(
          all_of(wdi_cov),
          ~ ifelse(
            is.na(.),
            round(mean(., na.rm = TRUE), 3),
            .
          )
        )
      ) %>%
      ungroup()
    
    # colnames(df_imputed)[colnames(df_imputed) == "country.y"] <- "country"
  }
  
  ##################################################################
  ############### prepare X ########################################
  ##################################################################
  
  #Get list of countries and resolutions
  if(include_hdi == 'yes'){
    countries <- sort(unique(df_imputed$country))
    
  }else{
    countries <- sort(unique(df_imputed$country.x))
    colnames(df_imputed)[colnames(df_imputed) == 'country.x'] <- 'country'
  }  

  n_country <- length(countries)
  
  if(include_hdi == 'yes'){
    cont_covariates = c(wdi_cov, "hdi","gdi")
  }else{
    cont_covariates = c(wdi_cov)
  }
  
  #Create country covariate matrix (n_country x 8)
  X_country <- matrix(0, nrow = n_country, ncol = length(cont_covariates))
  
  for(i in 1:n_country){
    c <- countries[i]
    
    # Subset rows for this country
    sub_df <- df_imputed[df_imputed$country == c, ] 
    
    ### first fill wdi covariates
    X_country[i, 1:length(wdi_cov)] <-
      sapply(wdi_cov, function(v) {
        x <- sub_df[[v]]
        ifelse(any(!is.na(x)), x[!is.na(x)][1], 0)
      })
    
    ## then fill hdi and gdi if include_hdi is "yes"
    if (include_hdi == 'yes'){
      X_country[i, (length(wdi_cov) + 1)] <- sub_df$hdi[!is.na(sub_df$hdi)][1]
      X_country[i, (length(wdi_cov) + 2)] <- sub_df$gdi[!is.na(sub_df$gdi)][1]
    }
  }   
  
  
  X <- X_country
  
  ## replace NA's by 0
  X[is.na(X)] <- 0
  
  # Assign column names
  if(include_hdi == 'yes'){
    colnames(X) <- c(wdi_cov, "hdi", "gdi")
  }else{
    colnames(X) <- c(wdi_cov)
  }
  
  # Assign row names to X
  rownames(X) <- c(countries)
  
  scaling <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  # scale only continuous covariates
  X[, cont_covariates] <- apply(
    X[, cont_covariates, drop = FALSE],
    2,
    scaling
  )
  ################ Adj prep ################################################
  
  df_a = df_imputed[,c('descr','country','vote')]
  
  df_a$vote_numeric <- ifelse(df_a$vote == "yes", 1.0,
                              ifelse(df_a$vote == "no", -1.0,
                                     ifelse(df_a$vote == "abstain", 0.0, NA)))
  df_a = subset(df_a, select = -vote)
  
  ##important!removing rows where descr ="" (year:1997)
  df_a <- df_a[df_a$descr != "", ]
  
  countries <- unique(df_a$country)
  descriptions <- c(unique(df_a$descr))
  
  ## important as 1997 description contain ""
  descriptions <- descriptions[descriptions != ""]
  
  
  # Initialize matrix with zeros
  vote_matrix <- matrix(0, nrow = length(countries), ncol = length(descriptions),
                        dimnames = list(countries, descriptions))
  
  # Fill matrix
  for (i in seq_len(nrow(df_a))) {
    ctry <- df_a$country[i]
    desc <- df_a$descr[i]
    
    vote_matrix[ctry, desc] <- df_a$vote_numeric[i]
  }
  
  
  ### !!!! VERY IMPORTANT TO MATCH THE ORDER ###
  # Reorder vote_matrix rows to match X
  vote_matrix <- vote_matrix[rownames(X), ]
  
  
  #### adj prep based on similarity measures###
  if (similarity_measure == 'cosine'){
    
    normalized_matrix = normalize(vote_matrix)
    
    Adj <- round(as.matrix(tcrossprod(normalized_matrix)),3)
    
  }else if (similarity_measure == 'jaccard'){
    
    Adj = outer(1:nrow(vote_matrix), 1:nrow(vote_matrix), Vectorize(function(i, j) jaccard_similarity(vote_matrix[i, ], vote_matrix[j, ])))
    
    Adj = round(as.matrix(Adj),3)
    
  }else {
    
    stop("Error!choose cosine or jaccard similarity")
  }
  
  ## no self-loop
  diag(Adj)<- 0
  Adj[Adj<0]<-0  # making the negative weights = 0 
  
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
  
  A = Adj_binary
  
  n = nrow(A)
  
  tau = mean(rowSums(A))
  
  inv_sqrt_deg <- ifelse(rowSums(A) > 0, 1 / sqrt(rowSums(A + tau * diag(n))), 0)
  
  D_inv_sqrt <- diag(inv_sqrt_deg)
  
  L_sym = tcrossprod(D_inv_sqrt %*% A, D_inv_sqrt) 
  
  L_sym[is.na(L_sym)] <- 0
  
  
  ## choose alpha ##
  LL    <- crossprod(L_sym)                       # L_tau L_tau  (= L^2, symmetric PSD)
  X_xt  <- tcrossprod(X); X_xt[is.na(X_xt)] <- 0   # XX^T
  
  lam1_LL <- eigs_sym(LL,   k = 1, which = "LM")$values[1]   # leading eigenvalue of L^2
  lam1_X  <- eigs_sym(X_xt, k = 1, which = "LM")$values[1]   # leading eigenvalue of XX^T
  
  alpha0  <- lam1_LL / lam1_X                      # Binkiewicz balancing value
  
  L_tilde = crossprod(L_sym) + alpha0*tcrossprod(X)

  eig <- RSpectra::eigs_sym(L_tilde,  k = K, which = "LM")  # "LM" = largest magnitude
  
  U <- eig$vectors  # n x K matrix
  
  row_norms <- sqrt(rowSums(U^2))
  row_norms[row_norms == 0] <- 1
  U_norm <- U / row_norms
  
  # Run k-means on rows of U_country
  set.seed(3456) 
  
  # Output
  kmeans_result <- kmeans(U_norm, centers = K, nstart = 10)
  
  kmeans_output_cluster <- kmeans_result$cluster
  
  country_clusters <- kmeans_output_cluster[1:n_country]
  
  # Map country names to clusters
  
  if (!identical(rownames(X), rownames(A))) {
    
    stop("X and A rownames mismatch.")
    
  }else{
    
    cluster_map = data.frame(
      country = rownames(A),          
      cluster = country_clusters,
      country_degrees = colSums(A))
  }
  
  
  return(cluster_map)
}


############################################################################
############## HOM ########################################################

## load helper function adj_prep_unga_r###########
source("adj_prep.R")


hom_unga_clustering <- function(input_year, similarity_measure, threshold = c('q1','q2','q3','mean'), K){
  
  
  A = adj_prep_from_unga_r(input_year = input_year, similarity_measure = 'cosine', threshold = threshold)
  
  n = nrow(A)
  
  tau = mean(rowSums(A))
  
  inv_sqrt_deg <- ifelse(rowSums(A) > 0, 1 / sqrt(rowSums(A + tau * diag(n))), 0)
  
  D_inv_sqrt <- diag(inv_sqrt_deg)
  
  L_sym = tcrossprod(D_inv_sqrt %*% A, D_inv_sqrt) 
  
  L_sym[is.na(L_sym)] <- 0
  
  eig <- RSpectra::eigs_sym(L_sym, k = K, which = "LM")  # "LM" = largest magnitude
  
  U <- eig$vectors  # n x K matrix
  
  row_norms <- sqrt(rowSums(U^2))
  row_norms[row_norms == 0] <- 1
  U_norm <- U / row_norms
  
  # Run k-means on rows of U_country
  set.seed(3456) #### reproducibility of cluster label 
  
  # Output
  kmeans_result <- kmeans(U_norm, centers = K, nstart = 10)
  
  # kmeans_result <- kmeans(U_norm, centers = K, nstart = 10)
  kmeans_output_cluster <- kmeans_result$cluster
  
  countries = rownames(A)
  country_degrees = colSums(A)
  
  ##################################################################
  ############ Plotting ############################################
  ##################################################################
  
  
  dfff = data.frame(country = c(countries),
                    country_degree = c(country_degrees),
                    cluster = kmeans_output_cluster,
                    size_of_nodes = c(0.1 * country_degrees + 5))
  
  dfff$country = toupper(dfff$country)
  
  
  if (!identical(toupper(row.names(dfff)), toupper(dfff$country))) {
    
    stop("Error: rownames(A) and countrynames are not identical in the same order.")
    
  }else{
    
    cluster_map <- dfff %>% dplyr::select(country, cluster, country_degree) ## mention dplyr::select-- because MASS also has select. There could be conflict
    
  }
  
  
  return(cluster_map)
  
}




