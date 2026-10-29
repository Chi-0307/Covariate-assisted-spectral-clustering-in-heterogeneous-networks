
##========DESCRIPTION OF THE CODE===============================================
## It loads the unga metric function. Unga metric function gives affiliation plots and LCP plots.
## It preloads the basic het-cov,hom-cov,het and hom functions, each of which produces a cluster map i.e a clustering of the countries into 3 clusters
## This code produces LCP plots for BRICS,SCO,LA,EU(for each of the methods).
## Also produces plots for Affiliation for the above mentioned organizations(for each of the methods).
## calculating affiliations and lc proportions for different years for different institutions(say, sco, eu, brics etc.)



library(future)
library(future.apply)
library(parallel)

## 
source("unga_metrics.R")
source('Helpers.R')


n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)


################# HOM-COV ######################################################

##### only gdp #################################################################

years <- 1990:2000

results_gdp_homcov_90_00 <- future_lapply(years, function(yr) {
  hom_cov_unga_final(input_year = yr, similarity_measure = 'cosine',
                     threshold = 'q2',
                     K = 3, 
                     wdi_indicators = c('NY.GDP.PCAP.CD'), 
                     include_hdi = 'no')
})

######### BRICS RESULTS ##########
brics_countries = toupper(c('BRAZIL', 'RUSSIA', 'INDIA', 'CHINA'))#,'SOUTH AFRICA'

years = c(1990:2000)

out_brics_homcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_brics_homcov_gdp = unga_proportion_metrics(member_country_list = brics_countries, 
                                                 cluster_map = results_gdp_homcov_90_00[[(yr-1990 + 1)]], 
                                                 affiliation_countries = c('RUSSIA','CHINA'))
  
  
  return(list(out_brics_homcov_gdp = out_brics_homcov_gdp))
  
})


cov_sets <- list(
  c("gdp")
)


### rus affiliation ###
rus_prop_list_homcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_homcov_", paste(covs, collapse = "_"))
  sapply(out_brics_homcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_homcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homcov")
)

### china affiliation ###
china_prop_list_homcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_homcov_", paste(covs, collapse = "_"))
  sapply(out_brics_homcov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_homcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("china_", paste(covs, collapse = "_"), "_homcov")
)

### largest cluster proportion #####
lc_prop_list_homcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_homcov_", paste(covs, collapse = "_"))
  sapply(out_brics_homcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_homcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homcov")
)


####### SCO #################################################################

sco_countries = c('CHINA', 'RUSSIA', 'KAZAKHSTAN', 'KYRGYZSTAN', 'TAJIKISTAN', 'UZBEKISTAN')

years = c(1990:2000)

out_sco_homcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_sco_homcov_gdp = unga_proportion_metrics(member_country_list = sco_countries, 
                                               cluster_map = results_gdp_homcov_90_00[[(yr-1990 + 1)]], 
                                               affiliation_countries = c('RUSSIA','CHINA'))
  
  
  return(list(out_sco_homcov_gdp = out_sco_homcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_homcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_homcov_", paste(covs, collapse = "_"))
  sapply(out_sco_homcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_homcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homcov")
)

### china affiliation ###
china_prop_list_homcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_homcov_", paste(covs, collapse = "_"))
  sapply(out_sco_homcov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_homcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("china_", paste(covs, collapse = "_"), "_homcov")
)

### largest cluster proportion #####
lc_prop_list_homcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_homcov_", paste(covs, collapse = "_"))
  sapply(out_sco_homcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_homcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homcov")
)

##########################################################################################
############# LA (latin american countries) ####################################

la_countries = c('ARGENTINA','BOLIVIA','BRAZIL','CHILE','COLOMBIA','ECUADOR','MEXICO',
                 'PARAGUAY','PERU','URUGUAY','VENEZUELA',
                 'CUBA',
                 'PANAMA')


years = c(1990:2000)

out_la_homcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_la_homcov_gdp = unga_proportion_metrics(member_country_list = la_countries, 
                                              cluster_map = results_gdp_homcov_90_00[[(yr-1990 + 1)]], 
                                              affiliation_countries = c('RUSSIA','CHINA'))
  
  
  return(list(out_la_homcov_gdp = out_la_homcov_gdp))
  
})



### largest cluster proportion #####
lc_prop_list_homcov_la <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_la_homcov_", paste(covs, collapse = "_"))
  sapply(out_la_homcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_homcov_la) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homcov")
)

############### EU ###########################################################

eu90_94_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK','UNITED STATES','Russia')
eu95_00_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK', 'Austria', 'Finland', 'Sweden','UNITED STATES','Russia')

######## 1990 - 1994 #######################
years_90_94 = 1990:1994

out_eu90_94_homcov = future_lapply(years_90_94, function(yr) {
  
  ### without hdi and gdi ######
  out_eu90_94_homcov_gdp = unga_proportion_metrics(member_country_list = eu90_94_countries, 
                                                   cluster_map = results_gdp_homcov_90_00[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  
  return(list(out_eu90_94_homcov_gdp = out_eu90_94_homcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_homcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_homcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_homcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homcov")
)

### UNITED STATES affiliation ###
usa_prop_list_homcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_homcov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_homcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("usa_", paste(covs, collapse = "_"), "_homcov")
)

### largest cluster proportion #####
lc_prop_list_homcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_homcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_homcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homcov")
)


######## 1995 - 2000 #######################
years_95_00 = 1995:2000

out_eu95_00_homcov = future_lapply(years_95_00, function(yr) {
  
  ### without hdi and gdi ######
  out_eu95_00_homcov_gdp = unga_proportion_metrics(member_country_list = eu95_00_countries, 
                                                   cluster_map = results_gdp_homcov_90_00[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  
  return(list(out_eu95_00_homcov_gdp = out_eu95_00_homcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_homcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_homcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_homcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homcov")
)

### UNITED STATES affiliation ###
usa_prop_list_homcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_homcov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_homcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("usa_", paste(covs, collapse = "_"), "_homcov")
)

### largest cluster proportion #####
lc_prop_list_homcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_homcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_homcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_homcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homcov")
)



# Stop parallelization
plan(sequential)

##############################################################
##################### HOM NO COV #############################

library(future)
library(future.apply)
library(parallel)


n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)



#### HOM-NO-COV ###################################

gdpout_nhnc <- future_lapply(years, function(yr){
  
  gdp <- hom_unga_clustering(input_year = yr, 
                             similarity_measure = 'cosine', 
                             threshold = c('q2'), 
                             K = 3)
  
  return(gdp)
  
}) 

############## sco ###################################
sco_countries = c('CHINA', 'RUSSIA', 'KAZAKHSTAN', 'KYRGYZSTAN', 'TAJIKISTAN', 'UZBEKISTAN')


out_sco_homnocov = future_lapply(years, function(yr) {
  
  
  out_sco_homnocov_gdp = unga_proportion_metrics(member_country_list = sco_countries, 
                                                 cluster_map = gdpout_nhnc[[(yr-1990 + 1)]], 
                                                 affiliation_countries = c('RUSSIA','CHINA'))
  
  return(list(out_sco_homnocov_gdp = out_sco_homnocov_gdp))
  
})



sets = c('gdp')

### rus affiliation ###
rus_prop_list_nhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_homnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_homnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_nhnc_sco) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homnocov")
)


### china affiliation ###
china_prop_list_nhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_homnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_homnocov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_nhnc_sco) <- sapply(
  sets,
  function(covs) paste0("ch_", paste(covs, collapse = "_"), "_homnocov")
)

### largest cluster proportion #####
lc_prop_list_nhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_homnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_homnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_nhnc_sco) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homnocov")
)


############## brics ###################################
brics_countries = toupper(c('BRAZIL', 'RUSSIA', 'INDIA', 'CHINA'))#,'SOUTH AFRICA'

out_brics_homnocov = future_lapply(years, function(yr) {
  
  
  out_brics_homnocov_gdp = unga_proportion_metrics(member_country_list = brics_countries, 
                                                   cluster_map = gdpout_nhnc[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','CHINA'))
  
  
  
  return(list(out_brics_homnocov_gdp = out_brics_homnocov_gdp))
  
})



sets = c('gdp')

### rus affiliation ###
rus_prop_list_nhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_homnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_homnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_nhnc_brics) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homnocov")
)


### china affiliation ###
china_prop_list_nhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_homnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_homnocov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_nhnc_brics) <- sapply(
  sets,
  function(covs) paste0("ch_", paste(covs, collapse = "_"), "_homnocov")
)


### largest cluster proportion #####
lc_prop_list_nhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_homnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_homnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_nhnc_brics) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homnocov")
)

############# LA ##############################################################

la_countries = c('ARGENTINA','BOLIVIA','BRAZIL','CHILE','COLOMBIA','ECUADOR','MEXICO',
                 'PARAGUAY','PERU','URUGUAY','VENEZUELA',
                 'CUBA',
                 'PANAMA')

out_la_homnocov = future_lapply(years, function(yr) {
  
  out_la_homnocov_gdp = unga_proportion_metrics(member_country_list = la_countries, 
                                                cluster_map = gdpout_nhnc[[(yr-1990 + 1)]], 
                                                affiliation_countries = c('RUSSIA','CHINA'))
  
  return(list(out_la_homnocov_gdp = out_la_homnocov_gdp))
})

sets = c('gdp')

### largest cluster proportion #####
lc_prop_list_nhnc_la <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_la_homnocov_", paste(covs, collapse = "_"))
  sapply(out_la_homnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_nhnc_la) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homnocov")
)

###################### eu #######################################################

eu90_94_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'United Kingdom','United States','Russia')
eu95_00_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'United Kingdom', 'Austria', 'Finland', 'Sweden','United States','Russia')

######## 1990 - 1994 #######################
years_90_94 = 1990:1994

out_eu_90_94_homnocov = future_lapply(years_90_94, function(yr) {
  
  
  out_eu_90_94_homnocov_gdp = unga_proportion_metrics(member_country_list = eu90_94_countries, 
                                                      cluster_map = gdpout_nhnc[[(yr-1990 + 1)]], 
                                                      affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  return(list(out_eu_90_94_homnocov_gdp = out_eu_90_94_homnocov_gdp))
  
})


sets = c('gdp')

### rus affiliation ###
rus_prop_list_nhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_homnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_nhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homnocov")
)


### usa affiliation ###
usa_prop_list_nhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_homnocov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_nhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("us_", paste(covs, collapse = "_"), "_homnocov")
)


### largest cluster proportion #####
lc_prop_list_nhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_homnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_nhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homnocov")
)

######## 1995 - 2000 #######################

years_95_00 = 1995:2000

out_eu_95_00_homnocov = future_lapply(years_95_00, function(yr) {
  
  
  ## covariate doesnt matter...so take any one of the combinations
  ## i.e result will be same accross all the cov combinations.
  out_eu_95_00_homnocov_gdp = unga_proportion_metrics(member_country_list = eu95_00_countries, 
                                                      cluster_map = gdpout_nhnc[[(yr-1995 + 1)]], 
                                                      affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  return(list(out_eu_95_00_homnocov_gdp = out_eu_95_00_homnocov_gdp))
  
})


sets = c('gdp')

### rus affiliation ###
rus_prop_list_nhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_homnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_nhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_homnocov")
)


### usa affiliation ###
usa_prop_list_nhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_homnocov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_nhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("us_", paste(covs, collapse = "_"), "_homnocov")
)


### largest cluster proportion #####
lc_prop_list_nhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_homnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_homnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_nhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_homnocov")
)

plan(sequential)
##################################################

rus_prop_list_nhnc_eu = c(c(rus_prop_list_nhnc_eu_90_94$rus_gdp_homnocov),
                          c(rus_prop_list_nhnc_eu_95_00$rus_gdp_homnocov))

usa_prop_list_nhnc_eu = c(usa_prop_list_nhnc_eu_90_94$us_gdp_homnocov,
                          usa_prop_list_nhnc_eu_95_00$us_gdp_homnocov)

lc_prop_list_nhnc_eu  = c(lc_prop_list_nhnc_eu_90_94$lc_gdp_homnocov,
                          lc_prop_list_nhnc_eu_95_00$lc_gdp_homnocov)

############################################################################
############################################################################
##################### HET - COV ############################################
############################################################################

library(future)
library(future.apply)
library(parallel)

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)

################################ HET-COV #######################################

##### only gdp #########################
years <- 1990:2000

results_gdp_hetcov_90_00 <- future_lapply(years, function(yr) {
  het_cov_unga_final(input_year = yr,
                     K = 3, 
                     method = c('het_yes_cov'), 
                     wdi_indicators = c('NY.GDP.PCAP.CD'), 
                     include_hdi = 'no')
})



################ BRICS RESULTS ##########
brics_countries = toupper(c('BRAZIL', 'RUSSIA', 'INDIA', 'CHINA'))#,'SOUTH AFRICA'

years = c(1990:2000)
out_brics_hetcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_brics_hetcov_gdp = unga_proportion_metrics(member_country_list = brics_countries, 
                                                 cluster_map = results_gdp_hetcov_90_00[[(yr-1990 + 1)]], 
                                                 affiliation_countries = c('RUSSIA','CHINA'))
  
  
  return(list(out_brics_hetcov_gdp = out_brics_hetcov_gdp))
  
})


cov_sets <- list(
  c("gdp")
)


### rus affiliation ###
rus_prop_list_hetcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_hetcov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_hetcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetcov")
)

### china affiliation ###
china_prop_list_hetcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_hetcov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetcov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_hetcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("china_", paste(covs, collapse = "_"), "_hetcov")
)

### largest cluster proportion #####
lc_prop_list_hetcov_brics <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_brics_hetcov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_hetcov_brics) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetcov")
)


####### SCO #################################################################

sco_countries = c('CHINA', 'RUSSIA', 'KAZAKHSTAN', 'KYRGYZSTAN', 'TAJIKISTAN', 'UZBEKISTAN')

years = c(1990:2000)

out_sco_hetcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_sco_hetcov_gdp = unga_proportion_metrics(member_country_list = sco_countries, 
                                               cluster_map = results_gdp_hetcov_90_00[[(yr-1990 + 1)]], 
                                               affiliation_countries = c('RUSSIA','CHINA'))
  
  
  
  return(list(out_sco_hetcov_gdp = out_sco_hetcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_hetcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_hetcov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_hetcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetcov")
)

### china affiliation ###
china_prop_list_hetcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_hetcov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetcov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_hetcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("china_", paste(covs, collapse = "_"), "_hetcov")
)

### largest cluster proportion #####
lc_prop_list_hetcov_sco <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_sco_hetcov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_hetcov_sco) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetcov")
)


############# LA (latin american countries) ####################################

la_countries = c('ARGENTINA','BOLIVIA','BRAZIL','CHILE','COLOMBIA','ECUADOR','MEXICO',
                 'PARAGUAY','PERU','URUGUAY','VENEZUELA',
                 'CUBA',
                 'PANAMA')


years = c(1990:2000)

out_la_hetcov = future_lapply(years, function(yr) {
  
  ### without hdi and gdi ######
  out_la_hetcov_gdp = unga_proportion_metrics(member_country_list = la_countries, 
                                              cluster_map = results_gdp_hetcov_90_00[[(yr-1990 + 1)]], 
                                              affiliation_countries = c('RUSSIA','CHINA'))## affiliation is redundant here
  
 
  return(list(out_la_hetcov_gdp = out_la_hetcov_gdp))
  
})



### largest cluster proportion #####
lc_prop_list_hetcov_la <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_la_hetcov_", paste(covs, collapse = "_"))
  sapply(out_la_hetcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_hetcov_la) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetcov")
)

############### EU ###########################################################

eu90_94_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK','USA','Russia')
eu95_00_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK', 'Austria', 'Finland', 'Sweden','USA','Russia')

######## 1990 - 1994 #######################
years_90_94 = 1990:1994

out_eu90_94_hetcov = future_lapply(years_90_94, function(yr) {
  
  ### without hdi and gdi ######
  out_eu90_94_hetcov_gdp = unga_proportion_metrics(member_country_list = eu90_94_countries, 
                                                   cluster_map = results_gdp_hetcov_90_00[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  
  return(list(out_eu90_94_hetcov_gdp = out_eu90_94_hetcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_hetcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_hetcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_hetcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetcov")
)

### USA affiliation ###
usa_prop_list_hetcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_hetcov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_hetcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("usa_", paste(covs, collapse = "_"), "_hetcov")
)

### largest cluster proportion #####
lc_prop_list_hetcov_eu90_94 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu90_94_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu90_94_hetcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_hetcov_eu90_94) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetcov")
)


######## 1995 - 2000 #######################
years_95_00 = 1995:2000

out_eu95_00_hetcov = future_lapply(years_95_00, function(yr) {
  
  ### without hdi and gdi ######
  out_eu95_00_hetcov_gdp = unga_proportion_metrics(member_country_list = eu95_00_countries, 
                                                   cluster_map = results_gdp_hetcov_90_00[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  
  return(list(out_eu95_00_hetcov_gdp = out_eu95_00_hetcov_gdp))
  
})


### rus affiliation ###
rus_prop_list_hetcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_hetcov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_hetcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetcov")
)

### USA affiliation ###
usa_prop_list_hetcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_hetcov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_hetcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("usa_", paste(covs, collapse = "_"), "_hetcov")
)

### largest cluster proportion #####
lc_prop_list_hetcov_eu95_00 <- future_lapply(cov_sets, function(covs) {
  obj_name <- paste0("out_eu95_00_hetcov_", paste(covs, collapse = "_"))
  sapply(out_eu95_00_hetcov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_hetcov_eu95_00) <- sapply(
  cov_sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetcov")
)
#####################################################################################################################
############################ HET ##############################
#### het-no-cov(yhnc) ###################################

## the covariates doesnt matter, taking any combination of covariates(not only gdp) would have given the same result)

gdpout_yhnc <- future_lapply(years, function(yr){
  
  gdp <- het_cov_unga_final(
    input_year = yr,
    K = 3,
    method = "het_no_cov",
    wdi_indicators = c("NY.GDP.PCAP.CD"),
    include_hdi = "no"
  )
  
  return(gdp)
  
}) 
####################################################################################################

############## sco ###################################
sco_countries = c('CHINA', 'RUSSIA', 'KAZAKHSTAN', 'KYRGYZSTAN', 'TAJIKISTAN', 'UZBEKISTAN')

out_sco_hetnocov = future_lapply(years, function(yr) {
  

  
  out_sco_hetnocov_gdp = unga_proportion_metrics(member_country_list = sco_countries, 
                                                 cluster_map = gdpout_yhnc[[(yr-1990 + 1)]], 
                                                 affiliation_countries = c('RUSSIA','CHINA'))
  
  
  return(list(out_sco_hetnocov_gdp = out_sco_hetnocov_gdp))
  
})



sets = c('gdp')

### rus affiliation ###
rus_prop_list_yhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_yhnc_sco) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetnocov")
)


### china affiliation ###
china_prop_list_yhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetnocov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_yhnc_sco) <- sapply(
  sets,
  function(covs) paste0("ch_", paste(covs, collapse = "_"), "_hetnocov")
)

### largest cluster proportion #####
lc_prop_list_yhnc_sco <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_sco_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_sco_hetnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_yhnc_sco) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetnocov")
)
#### check! the results doesnt depend on covariates.

############## brics ###################################
brics_countries = toupper(c('BRAZIL', 'RUSSIA', 'INDIA', 'CHINA'))#,'SOUTH AFRICA'

out_brics_hetnocov = future_lapply(years, function(yr) {
  
  
  
  out_brics_hetnocov_gdp = unga_proportion_metrics(member_country_list = brics_countries, 
                                                   cluster_map = gdpout_yhnc[[(yr-1990 + 1)]], 
                                                   affiliation_countries = c('RUSSIA','CHINA'))
  
  
  
  return(list(out_brics_hetnocov_gdp = out_brics_hetnocov_gdp))
  
})



sets = c('gdp')

### rus affiliation ###
rus_prop_list_yhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_yhnc_brics) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetnocov")
)


### china affiliation ###
china_prop_list_yhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetnocov, function(x)
    x[[obj_name]]$prop_ch
  )
})

names(china_prop_list_yhnc_brics) <- sapply(
  sets,
  function(covs) paste0("ch_", paste(covs, collapse = "_"), "_hetnocov")
)

### largest cluster proportion #####
lc_prop_list_yhnc_brics <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_brics_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_brics_hetnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_yhnc_brics) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetnocov")
)

############# LA ##############################################################


la_countries = c('ARGENTINA','BOLIVIA','BRAZIL','CHILE','COLOMBIA','ECUADOR','MEXICO',
                 'PARAGUAY','PERU','URUGUAY','VENEZUELA',
                 'CUBA',
                 'PANAMA')

out_la_hetnocov = future_lapply(years, function(yr) {
  
  ## covariate doesnt matter...so take any one of the combinations
  ## i.e result will be same accross all the cov combinations.
  out_la_hetnocov_gdp = unga_proportion_metrics(member_country_list = la_countries, 
                                                    cluster_map = gdpout_yhnc[[(yr-1990 + 1)]], 
                                                    affiliation_countries = c('RUSSIA','CHINA'))
  
  return(list(out_la_hetnocov_gdp = out_la_hetnocov_gdp))
})

sets = c('gdp')

### largest cluster proportion #####
lc_prop_list_yhnc_la <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_la_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_la_hetnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_yhnc_la) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetnocov")
)

###################### eu #######################################################

eu90_94_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK','UNITED STATES','Russia')
eu95_00_countries = c('Belgium', 'Denmark', 'France', 'Germany', 'Greece', 'Ireland', 'Italy', 'Luxembourg', 'Netherlands', 'Portugal', 'Spain', 'UK', 'Austria', 'Finland', 'Sweden','UNITED STATES','Russia')

######## 1990 - 1994 #######################
years_90_94 = 1990:1994

out_eu_90_94_hetnocov = future_lapply(years_90_94, function(yr) {
  
  ## result will be same accross all the cov combinations.
  out_eu_90_94_hetnocov_gdp = unga_proportion_metrics(member_country_list = eu90_94_countries, 
                                                          cluster_map = gdpout_yhnc[[(yr-1990 + 1)]], 
                                                          affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  return(list(out_eu_90_94_hetnocov_gdp = out_eu_90_94_hetnocov_gdp))
  
})


sets = c('gdp')

### rus affiliation ###
rus_prop_list_yhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_hetnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_yhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetnocov")
)


### UNITED STATES affiliation ###
usa_prop_list_yhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_hetnocov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_yhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("us_", paste(covs, collapse = "_"), "_hetnocov")
)


### largest cluster proportion #####
lc_prop_list_yhnc_eu_90_94 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_90_94_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_90_94_hetnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_yhnc_eu_90_94) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetnocov")
)


######## 1995 - 2000 #######################
years_95_00 = 1995:2000

out_eu_95_00_hetnocov = future_lapply(years_95_00, function(yr) {
  
  
  ## result will be same accross all the cov combinations.
  out_eu_95_00_hetnocov_gdp = unga_proportion_metrics(member_country_list = eu95_00_countries, 
                                                          cluster_map = gdpout_yhnc[[(yr-1995 + 1)]], 
                                                          affiliation_countries = c('RUSSIA','UNITED STATES'))
  
  return(list(out_eu_95_00_hetnocov_gdp = out_eu_95_00_hetnocov_gdp))
  
})


sets = c('gdp')

### rus affiliation ###
rus_prop_list_yhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_hetnocov, function(x)
    x[[obj_name]]$prop_rus
  )
})

names(rus_prop_list_yhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("rus_", paste(covs, collapse = "_"), "_hetnocov")
)


### usa affiliation ###
usa_prop_list_yhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_hetnocov, function(x)
    x[[obj_name]]$prop_us
  )
})

names(usa_prop_list_yhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("us_", paste(covs, collapse = "_"), "_hetnocov")
)


### largest cluster proportion #####
lc_prop_list_yhnc_eu_95_00 <- future_lapply(sets, function(covs) {
  obj_name <- paste0("out_eu_95_00_hetnocov_", paste(covs, collapse = "_"))
  sapply(out_eu_95_00_hetnocov, function(x)
    x[[obj_name]]$prop_lc
  )
})

names(lc_prop_list_yhnc_eu_95_00) <- sapply(
  sets,
  function(covs) paste0("lc_", paste(covs, collapse = "_"), "_hetnocov")
)

##################################################

rus_prop_list_yhnc_eu = c(c(rus_prop_list_yhnc_eu_90_94$rus_gdp_hetnocov),
                          c(rus_prop_list_yhnc_eu_95_00$rus_gdp_hetnocov))

usa_prop_list_yhnc_eu = c(usa_prop_list_yhnc_eu_90_94$us_gdp_hetnocov,
                          usa_prop_list_yhnc_eu_95_00$us_gdp_hetnocov)

lc_prop_list_yhnc_eu  = c(lc_prop_list_yhnc_eu_90_94$lc_gdp_hetnocov,
                          lc_prop_list_yhnc_eu_95_00$lc_gdp_hetnocov)
plan(sequential)

################################################################################
######################## PLOTS #################################################
################################################################################

years <- 1990:2000
path = "Plots/lcp_and_aff_plots"

la_path = path;sco_path=path;brics_path=path;eu_path=path
################ LA PLOTS ######################################################

######### lc ###################################################################
largest_cl_yhnc_la = c(lc_prop_list_yhnc_la$lc_gdp_hetnocov) 

largest_cl_nhnc_la = c(lc_prop_list_nhnc_la$lc_gdp_homnocov)

future_lapply(names(lc_prop_list_hetcov_la), function(nm) {
  
  pdf(
    file = file.path(la_path, paste0("largest_cluster_", nm, "_la.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)

  
  plot(
    years, lc_prop_list_hetcov_la[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Largest Cluster Proportion",
    main = ""
  )
  
  lines(years, lc_prop_list_homcov_la[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, largest_cl_yhnc_la,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, largest_cl_nhnc_la,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})


################ BRICS PLOTS ######################################################

## largest cluster ####

largest_cl_yhnc_brics = c(lc_prop_list_yhnc_brics$lc_gdp_hetnocov)#c(3/4,3/4,3/4,3/4,3/4,3/4,2/4,2/4,2/4,2/4,2/4)
largest_cl_nhnc_brics = c(lc_prop_list_nhnc_brics$lc_gdp_homnocov)#c(4/4,3/4,2/4,3/4,3/4,2/4,3/4,2/4,3/4,3/4,3/4)

future_lapply(names(lc_prop_list_hetcov_brics), function(nm) {
  
  pdf(
    file = file.path(brics_path, paste0("largest_cluster_", nm, "_brics.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)

  
  plot(
    years, lc_prop_list_hetcov_brics[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Largest Cluster Proportion",
    main = ""
  )
  
  lines(years, lc_prop_list_homcov_brics[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, largest_cl_yhnc_brics,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, largest_cl_nhnc_brics,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

#### ch_aff #####################


ch_aff_yhnc_brics = c(china_prop_list_yhnc_brics$ch_gdp_hetnocov) ## choose any one 
ch_aff_nhnc_brics = c(china_prop_list_nhnc_brics$ch_gdp_homnocov)

future_lapply(names(china_prop_list_hetcov_brics), function(nm) {
  
  pdf(
    file = file.path(brics_path, paste0("ch_aff_", nm, "_brics.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)

  
  plot(
    years, china_prop_list_hetcov_brics[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(China)",
    main = " "
  )
  
  lines(years, china_prop_list_homcov_brics[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, china_prop_list_yhnc_brics$ch_gdp_hetnocov,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, china_prop_list_nhnc_brics$ch_gdp_homnocov,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

###### rus_aff ####################


rus_aff_yhnc_brics = c(rus_prop_list_yhnc_brics$rus_gdp_hetnocov) ## choose any one 
rus_aff_nhnc_brics = c(rus_prop_list_nhnc_brics$rus_gdp_homnocov)

future_lapply(names(rus_prop_list_hetcov_brics), function(nm) {
  
  pdf(
    file = file.path(brics_path, paste0("rus_aff_", nm, "_brics.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)

  
  plot(
    years, rus_prop_list_hetcov_brics[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(Russia)",
    main = " "
  )
  
  lines(years, rus_prop_list_homcov_brics[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, rus_prop_list_yhnc_brics$rus_gdp_hetnocov,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, rus_prop_list_nhnc_brics$rus_gdp_homnocov,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})



################################ SCO Plots #############################################


##### lc ##################################################

largest_cl_yhnc_sco = c(lc_prop_list_yhnc_sco$lc_gdp_hetnocov) 
largest_cl_nhnc_sco = c(lc_prop_list_nhnc_sco$lc_gdp_homnocov) 


future_lapply(names(lc_prop_list_hetcov_sco), function(nm) {
  
  pdf(
    file = file.path(sco_path, paste0("largest_cluster", nm, "_sco.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- range(
    lc_prop_list_hetcov_sco[[nm]],
    lc_prop_list_hetcov_sco[[gsub("hetcov", "homcov", nm)]],
    largest_cl_yhnc_sco,
    largest_cl_nhnc_sco
  )
  
  plot(
    years, lc_prop_list_hetcov_sco[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Largest Cluster Proportion",
    main = ""
  )
  
  lines(years, lc_prop_list_homcov_sco[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, largest_cl_yhnc_sco,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, largest_cl_nhnc_sco,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

#### ch_aff #####################

ch_aff_yhnc_sco = c(china_prop_list_yhnc_sco$ch_gdp_hetnocov) ## choose any one 
ch_aff_nhnc_sco = c(china_prop_list_nhnc_sco$ch_gdp_homnocov)

future_lapply(names(china_prop_list_hetcov_sco), function(nm) {
  
  pdf(
    file = file.path(sco_path, paste0("ch_aff_", nm, "_sco.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)
  
  plot(
    years, china_prop_list_hetcov_sco[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(China)",
    main = ""
  )
  
  lines(years, china_prop_list_homcov_sco[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, china_prop_list_yhnc_sco$ch_gdp_hetnocov,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, china_prop_list_nhnc_sco$ch_gdp_homnocov,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

###### rus_aff ####################

rus_aff_yhnc_sco = c(rus_prop_list_yhnc_sco$rus_gdp_hetnocov) 
rus_aff_nhnc_sco = c(rus_prop_list_nhnc_sco$rus_gdp_homnocov)

future_lapply(names(rus_prop_list_hetcov_sco), function(nm) {
  
  pdf(
    file = file.path(sco_path, paste0("rus_aff_", nm, "_sco.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)
  
  plot(
    years, rus_prop_list_hetcov_sco[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(Russia)",
    main = ""
  )
  
  lines(years, rus_prop_list_homcov_sco[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, rus_prop_list_yhnc_sco$rus_gdp_hetnocov,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, rus_prop_list_nhnc_sco$rus_gdp_homnocov,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})


################## EU Plots #############################################

##### lc ##################################################
largest_cl_yhnc_eu = c(lc_prop_list_yhnc_eu) 
largest_cl_nhnc_eu = c(lc_prop_list_nhnc_eu) 

### joining eu-90-94 and eu-95-00 results for hetcov 
names <- union(names(lc_prop_list_hetcov_eu90_94),
               names(lc_prop_list_hetcov_eu95_00))

lc_prop_list_hetcov_eu90_00 <- setNames(
  lapply(names, \(n) c(lc_prop_list_hetcov_eu90_94[[n]],
                       lc_prop_list_hetcov_eu95_00[[n]])),
  names
)

### joining eu-90-94 and eu-95-00 results for homcov 
names <- union(names(lc_prop_list_homcov_eu90_94),
               names(lc_prop_list_homcov_eu95_00))

lc_prop_list_homcov_eu90_00 <- setNames(
  lapply(names, \(n) c(lc_prop_list_homcov_eu90_94[[n]],
                       lc_prop_list_homcov_eu95_00[[n]])),
  names
)

future_lapply(names(lc_prop_list_hetcov_eu90_00), function(nm) {
  
  pdf(
    file = file.path(eu_path, paste0("largest_cluster", nm, "_eu.pdf")),
    width = 7, height = 5
  )
  
  
  ylim_all <- c(0,1)
  
  plot(
    years, lc_prop_list_hetcov_eu90_00[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Largest Cluster Proportion",
    main = ""
  )
  
  lines(years, lc_prop_list_homcov_eu90_00[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, largest_cl_yhnc_eu,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, largest_cl_nhnc_eu,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

#### us_aff #####################

usa_prop_list_yhnc_eu = c(usa_prop_list_yhnc_eu)
usa_prop_list_nhnc_eu = c(usa_prop_list_nhnc_eu)

### joining eu-90-94 and eu-95-00 results for hetcov 
names <- union(names(usa_prop_list_hetcov_eu90_94),
               names(usa_prop_list_hetcov_eu95_00))

usa_prop_list_hetcov_eu90_00 <- setNames(
  lapply(names, \(n) c(usa_prop_list_hetcov_eu90_94[[n]],
                       usa_prop_list_hetcov_eu95_00[[n]])),
  names
)

### joining eu-90-94 and eu-95-00 results for homcov 
names <- union(names(usa_prop_list_homcov_eu90_94),
               names(usa_prop_list_homcov_eu95_00))

usa_prop_list_homcov_eu90_00 <- setNames(
  lapply(names, \(n) c(usa_prop_list_homcov_eu90_94[[n]],
                       usa_prop_list_homcov_eu95_00[[n]])),
  names
)



future_lapply(names(usa_prop_list_hetcov_eu90_00), function(nm) {
  
  pdf(
    file = file.path(eu_path, paste0("us_aff_", nm, "_eu.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)
  
  plot(
    years, usa_prop_list_hetcov_eu90_00[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(usa)",
    main = ""
  )
  
  lines(years, usa_prop_list_homcov_eu90_00[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, usa_prop_list_yhnc_eu,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, usa_prop_list_nhnc_eu,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})

###### rus_aff ####################

rus_prop_list_yhnc_eu = c(rus_prop_list_yhnc_eu)
rus_prop_list_nhnc_eu = c(rus_prop_list_nhnc_eu)

### joining eu-90-94 and eu-95-00 results for hetcov 
names <- union(names(rus_prop_list_hetcov_eu90_94),
               names(rus_prop_list_hetcov_eu95_00))

rus_prop_list_hetcov_eu90_00 <- setNames(
  lapply(names, \(n) c(rus_prop_list_hetcov_eu90_94[[n]],
                       rus_prop_list_hetcov_eu95_00[[n]])),
  names
)

### joining eu-90-94 and eu-95-00 results for homcov 
names <- union(names(rus_prop_list_homcov_eu90_94),
               names(rus_prop_list_homcov_eu95_00))

rus_prop_list_homcov_eu90_00 <- setNames(
  lapply(names, \(n) c(rus_prop_list_homcov_eu90_94[[n]],
                       rus_prop_list_homcov_eu95_00[[n]])),
  names
)



future_lapply(names(rus_prop_list_hetcov_eu90_00), function(nm) {
  
  pdf(
    file = file.path(eu_path, paste0("rus_aff_", nm, "_eu.pdf")),
    width = 7, height = 5
  )
  
  ylim_all <- c(0,1)
  
  plot(
    years, rus_prop_list_hetcov_eu90_00[[nm]],
    type = "o",
    pch = 16,
    lwd = 2,
    col = "red3",
    ylim = ylim_all,
    xlab = "Year",
    ylab = "Affiliation(rus)",
    main = ""
  )
  
  lines(years, rus_prop_list_homcov_eu90_00[[gsub("hetcov", "homcov", nm)]],
        type = "o", pch = 7, lwd = 2, col = "royalblue1")
  
  lines(years, rus_prop_list_yhnc_eu,
        type = "o", pch = 15, lwd = 2, col = "tan4")
  
  lines(years, rus_prop_list_nhnc_eu,
        type = "o", pch = 18, lwd = 2, col = "limegreen")
  
  legend(
    "bottomleft",
    legend = c("Het-Cov", "Hom-Cov", "Het", "Hom"),
    col = c("red3", "royalblue1", "tan4", "limegreen"),
    pch = c(16, 17, 15, 18),
    lwd = 1.2,
    cex = 0.7,
    x.intersp = 0.5,
    y.intersp = 0.7,
    bty = "n"
  )
  
  dev.off()
  
  invisible(NULL)
})


# Stop parallelization
plan(sequential)

