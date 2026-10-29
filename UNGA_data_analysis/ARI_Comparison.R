
### Description #####
#this code compares the ari between the clusterings of two consecutive years, for each method.
#yhyc-yeshetyescov i.e  hetcov, yhnc-yeshetnocov i.e. het_only
#nhyc-nohetyescov  i.e. homcov, nhnc-nohetnocov  i.e  hom_only

library(unvotes)


wdi_ind_list = c("NY.GDP.PCAP.CD")


source('Helpers.R')

ari_comparison_between_years <-function(year1,year2,K){
  
  
  wdi_ind_list = c("NY.GDP.PCAP.CD")
  
  ## yhyc ##
  cl_yhyc_yr1 = het_cov_unga_final(input_year = year1, K=K, 
                                   method = 'het_yes_cov', wdi_indicators = wdi_ind_list, include_hdi='no')
  cl_yhyc_yr2 = het_cov_unga_final(input_year = year2, K=K, 
                                   method = 'het_yes_cov', wdi_indicators = wdi_ind_list, include_hdi='no')
  
  common_countries_yhyc = intersect(cl_yhyc_yr1$country,cl_yhyc_yr2$country)
  cl1_yhyc_common <- cl_yhyc_yr1[cl_yhyc_yr1$country %in% common_countries_yhyc, ]
  cl2_yhyc_common <- cl_yhyc_yr2[cl_yhyc_yr2$country %in% common_countries_yhyc, ]
  
  cl1_yhyc_common <- cl1_yhyc_common[order(cl1_yhyc_common$country), ]
  cl2_yhyc_common <- cl2_yhyc_common[order(cl2_yhyc_common$country), ]
  
  ari_yhyc <- adjustedRandIndex(cl1_yhyc_common$cluster, cl2_yhyc_common$cluster)  
  
  
  ### yhnc ###
  cl_yhnc_yr1 = het_cov_unga_final(input_year = year1, K=K, 
                                   method = 'het_no_cov',wdi_indicators=wdi_ind_list, include_hdi='no')
  cl_yhnc_yr2 = het_cov_unga_final(input_year = year2, K=K, 
                                   method = 'het_no_cov',wdi_indicators=wdi_ind_list, include_hdi='no')
  
  common_countries_yhnc = intersect(cl_yhnc_yr1$country,cl_yhnc_yr2$country)
  cl1_yhnc_common <- cl_yhnc_yr1[cl_yhnc_yr1$country %in% common_countries_yhnc, ]
  cl2_yhnc_common <- cl_yhnc_yr2[cl_yhnc_yr2$country %in% common_countries_yhnc, ]
  
  cl1_yhnc_common <- cl1_yhnc_common[order(cl1_yhnc_common$country), ]
  cl2_yhnc_common <- cl2_yhnc_common[order(cl2_yhnc_common$country), ]
  
  ari_yhnc <- adjustedRandIndex(cl1_yhnc_common$cluster, cl2_yhnc_common$cluster)  
  
  
  ## nhyc ##
  cl_nhyc_yr1 = hom_cov_unga_final(input_year = year1, 
                                   similarity_measure='cosine', threshold = c('q2'), K=K, wdi_indicators=wdi_ind_list, include_hdi='no')
  cl_nhyc_yr2 = hom_cov_unga_final(input_year = year2, 
                                   similarity_measure='cosine', threshold = c('q2'), K=K, wdi_indicators=wdi_ind_list, include_hdi='no')
  common_countries_nhyc = intersect(cl_nhyc_yr1$country,cl_nhyc_yr2$country)
  cl1_nhyc_common <- cl_nhyc_yr1[cl_nhyc_yr1$country %in% common_countries_nhyc, ]
  cl2_nhyc_common <- cl_nhyc_yr2[cl_nhyc_yr2$country %in% common_countries_nhyc, ]
  
  cl1_nhyc_common <- cl1_nhyc_common[order(cl1_nhyc_common$country), ]
  cl2_nhyc_common <- cl2_nhyc_common[order(cl2_nhyc_common$country), ]
  
  ari_nhyc <- adjustedRandIndex(cl1_nhyc_common$cluster, cl2_nhyc_common$cluster)  
  

  ## nhnc ##
  cl_nhnc_yr1 = hom_unga_clustering(input_year = year1, 
                                    similarity_measure = 'cosine', threshold = c('q2'), K = K)
  cl_nhnc_yr2 = hom_unga_clustering(input_year = year2, 
                                    similarity_measure = 'cosine', threshold = c('q2'), K = K)
  
  common_countries_nhnc = intersect(cl_nhnc_yr1$country,cl_nhnc_yr2$country)
  cl1_nhnc_common <- cl_nhnc_yr1[cl_nhnc_yr1$country %in% common_countries_nhnc, ]
  cl2_nhnc_common <- cl_nhnc_yr2[cl_nhnc_yr2$country %in% common_countries_nhnc, ]
  
  cl1_nhnc_common <- cl1_nhnc_common[order(cl1_nhnc_common$country), ]
  cl2_nhnc_common <- cl2_nhnc_common[order(cl2_nhnc_common$country), ]
  
  ari_nhnc <- adjustedRandIndex(cl1_nhnc_common$cluster, cl2_nhnc_common$cluster)  
  

  return(list(ari_yhyc = ari_yhyc,ari_yhnc = ari_yhnc,
              ari_nhyc = ari_nhyc,ari_nhnc = ari_nhnc))

}


################################################################################
library(future)
library(future.apply)
library(parallel)

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)


years <- 1990:1999

out_ari <- future_lapply(years, function(yr) {
    ari_comparison_between_years(yr, yr + 1, K = 3)
})

names(out_ari) <- paste0(years, "-", years + 1)
  
plan(sequential)


ari_mat <- do.call(cbind, out_ari)

# Set row names as year pairs
colnames(ari_mat) <- paste0(years, "-", years + 1)

# Set column names as methods
rownames(ari_mat) <- c("het-cov", "het", "hom-cov", "hom")

ari_df <- as.data.frame(ari_mat)
ari_df


