

###### affiliation and largest cluster proportion #########################################################


unga_proportion_metrics <- function(member_country_list, cluster_map, affiliation_countries){
  
  out1 <- cluster_map[toupper(cluster_map$country) %in% toupper(member_country_list), ]
  out1$country <- toupper(out1$country)
  
  split = split(out1$country,out1$cluster)
  
  ## largest cluster prop ##
  ## member countries ####
  members <- toupper(member_country_list)
  
  ## Step 1: identify largest cluster by size
  cluster_sizes <- sapply(split, length)
  largest_cluster <- split[[which.max(cluster_sizes)]]
  
  ## Step 2:  proportion in the largest cluster
  prop_largest_cluster <- sum(largest_cluster %in% members) / length(members)
  
  
  ## affiliation ###
  ## Identify cluster containing affiliationg countries (china, russia, usa etc.)
  affiliation_countries<-toupper(affiliation_countries)
  
  
  ## Identify cluster containing USA, China and Russia
  idx_us <- which(sapply(split, function(cl) "UNITED STATES" %in% cl)) 
  #### !!! mind that you should not write 'USA' instead of 'UNITED STATES' anywhere in the code

  idx_china <- which(sapply(split, function(cl) "CHINA" %in% cl))
  
  idx_rus   <- which(sapply(split, function(cl) "RUSSIA" %in% cl))
  
  cluster_us <- if (length(idx_us) > 0) split[[idx_us]] else NULL
  cluster_china <- if (length(idx_china) > 0) split[[idx_china]] else NULL
  cluster_rus  <- if (length(idx_rus) > 0)   split[[idx_rus]]   else NULL
  
  ## Proportion with usa
  if ("UNITED STATES" %in% affiliation_countries == TRUE){  ### dont write usa!
    prop_usa <- if (!is.null(cluster_us)) {
      sum(members %in% cluster_us) / length(members)
    } else 0
  }else{
    prop_usa <-"not_applicable"
  }

  
  ## Proportion with china
  if ("CHINA" %in% affiliation_countries == TRUE){
    prop_china <- if (!is.null(cluster_china)) {
      sum(members %in% cluster_china) / length(members)
    } else 0
  }else{
    prop_china <-"not_applicable"
  }
  
  
  ## Proportion with Russia
  if ("RUSSIA" %in% affiliation_countries == TRUE){
    prop_russia <- if (!is.null(cluster_rus)) {
      sum(members %in% cluster_rus) / length(members)
    } else 0
  }else{
    prop_russia <-"not_applicable"
  }

  return(list(prop_us = prop_usa, prop_ch = prop_china, prop_rus = prop_russia, prop_lc = prop_largest_cluster))  

}

 