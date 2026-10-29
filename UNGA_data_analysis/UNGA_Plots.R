################################################################################
### Description: ###############################################################
### This code give cluster plot and world clusters upon taking argument of cluster map
### 1990-2000 cluster plots are geenrated
################################################################################

source("Helpers.R")


unga_plots <- function(cluster_map, input_year){
  
  ##################################################################
  ######################## DATA PREP ###############################
  ##################################################################
  
  dfff = data.frame(
    country = c(cluster_map$country),
    country_degree = c(cluster_map$country_degrees),
    cluster = c(cluster_map$cluster),
    size_of_nodes = c(0.1 * c(cluster_map$country_degrees) + 5)
  )
  
  cluster_ids <- sort(unique(dfff$cluster))
  K <- length(cluster_ids)
  
  cluster_positions <- data.frame(
    cluster = cluster_ids,
    base_x = seq(0, 5 * (K - 1), by = 5),
    base_y = rep(0, K)
  )
  
  dfff <- dfff %>%
    left_join(cluster_positions, by = "cluster") %>%
    mutate(
      x = base_x + rnorm(n(), 0, 0.6),
      y = base_y + rnorm(n(), 0, 0.6)
    )
  
  ##################################################################
  ######################## COLOR SYSTEM ############################
  ##################################################################
  
  # Extract China and USA clusters
  usa_cluster   <- dfff$cluster[dfff$country == "United States"][1]
  china_cluster <- dfff$cluster[dfff$country == "China"][1]
  
  if (length(usa_cluster) == 0) usa_cluster <- NA
  if (length(china_cluster) == 0) china_cluster <- NA
  
  cluster_colors <- c()
  legend_labels  <- c()
  
  ##################################################################
  # CASE 1: China and USA in same cluster
  ##################################################################
  
  if (!is.na(china_cluster) && 
      !is.na(usa_cluster) && 
      china_cluster == usa_cluster) {
    
    message("china usa in same cluster")
    
    cluster_colors[as.character(china_cluster)] <- "lightskyblue"
    legend_labels[as.character(china_cluster)]  <- "Cluster A (China & USA)"
    
    special_clusters <- china_cluster
    
  } else {
    
    ################################################################
    # CASE 2: Separate clusters
    ################################################################
    
    special_clusters <- c()
    
    if (!is.na(china_cluster)) {
      cluster_colors[as.character(china_cluster)] <- "firebrick1"
      legend_labels[as.character(china_cluster)]  <- "Cluster A (China)"
      special_clusters <- c(special_clusters, china_cluster)
    }
    
    if (!is.na(usa_cluster)) {
      cluster_colors[as.character(usa_cluster)] <- "royalblue1"
      legend_labels[as.character(usa_cluster)]  <- "Cluster B (USA)"
      special_clusters <- c(special_clusters, usa_cluster)
    }
  }
  
  ##################################################################
  # Among remaining clusters - largest = green
  ##################################################################
  
  cluster_sizes <- sort(table(dfff$cluster), decreasing = TRUE)
  
  remaining_ranked <- cluster_sizes[
    !names(cluster_sizes) %in% as.character(special_clusters)
  ]
  
  largest_remaining_cluster <- NA
  
  if (length(remaining_ranked) >= 1) {
    largest_remaining_cluster <- as.numeric(names(remaining_ranked)[1])
    
    cluster_colors[as.character(largest_remaining_cluster)] <- "lawngreen"
    legend_labels[as.character(largest_remaining_cluster)]  <- "Cluster C"
  }
  
  ##################################################################
  # Remaining clusters → D, E, F ...
  ##################################################################
  
  other_clusters <- setdiff(cluster_ids,
                            c(special_clusters,
                              largest_remaining_cluster))
  
  extra_palette <- c(
    "gold1","magenta1","darkorange1",
    "turquoise1","purple1",
    "deeppink1","chartreuse1",
    "cyan1","hotpink1"
  )
  
  if (length(other_clusters) >= 1) {
    for (i in seq_along(other_clusters)) {
      cluster_colors[as.character(other_clusters[i])] <- extra_palette[i]
      legend_labels[as.character(other_clusters[i])]  <- 
        paste0("Cluster ", LETTERS[i + 3])
    }
  }
  
  ##################################################################
  # Legend Order: Red/Lightblue - Blue - Green - Others
  ##################################################################
  
  ordered_clusters <- c(
    special_clusters,
    largest_remaining_cluster,
    other_clusters
  )
  
  ordered_clusters <- ordered_clusters[!is.na(ordered_clusters)]
  
  cluster_colors <- cluster_colors[as.character(ordered_clusters)]
  legend_labels  <- legend_labels[as.character(ordered_clusters)]
  
  legend_levels <- ordered_clusters
  
  ########################################################
  ############## CLUSTER PLOT ############################
  ########################################################
  
  clustering_plot = ggplot(
    dfff,
    aes(x = x, y = y)
  ) +
    geom_point(
      aes(color = factor(cluster, levels = legend_levels),
          size = size_of_nodes),  # <-- controls size
      alpha = 0.8
    ) +
    geom_text_repel(
      aes(label = ifelse(country %in%
                           c("United States","United Kingdom","India",
                             "China","Israel","Palestine","Russia",
                             "Pakistan","Iran"), NA, country)),
      size = 2
    ) +
    geom_text_repel(
      data = subset(dfff,
                    country %in%
                      c("United States","United Kingdom","India",
                        "China","Israel","Palestine","Russia",
                        "Pakistan","Iran")),
      aes(label = country),
      fontface = "bold",
      size = 3
    ) +
    scale_color_manual(
      values = cluster_colors,
      labels = legend_labels,
      name = "Clusters"
    ) +
    scale_size_identity() +
    theme_void() +
    labs(
      title = paste0("Clusters: ", input_year),
    )
  
  ##################################################################
  ######################## WORLD MAP ###############################
  ##################################################################
  
  world <- map_data("world")
  
  # KEEPING ALL YOUR IMPORTANT NAME FIXES
  cluster_map$country[cluster_map$country == "United States"] <- "USA"
  cluster_map$country[cluster_map$country == "United Kingdom"] <- "UK"
  cluster_map$country[cluster_map$country ==  "Bosnia & Herzegovina"] <- "Bosnia and Herzegovina"
  cluster_map$country[cluster_map$country ==  "Czechia"] <- "Czech Republic"
  cluster_map$country[cluster_map$country ==  "Myanmar (Burma)"] <- "Myanmar"
  cluster_map$country[cluster_map$country ==  "Micronesia (Federated States of)"] <- "Micronesia"
  cluster_map$country[cluster_map$country ==  "São Tomé & Príncipe"] <- "Sao Tome and Principe"
  cluster_map$country[cluster_map$country ==  "St. Kitts & Nevis"] <- "Saint Kitts"
  cluster_map$country[cluster_map$country ==  "Antigua & Barbuda"] <- "Antigua"
  cluster_map$country[cluster_map$country ==  "Congo - Brazzaville"] <- "Democratic Republic of the Congo"
  cluster_map$country[cluster_map$country ==  "St. Vincent & Grenadines"] <- "Saint Vincent"
  cluster_map$country[cluster_map$country ==  "Eswatini"] <- "Swaziland"
  cluster_map$country[cluster_map$country ==  "St. Lucia"] <- "Saint Lucia"
  cluster_map$country[cluster_map$country ==  "Trinidad & Tobago"] <- "Trinidad"
  cluster_map$country[cluster_map$country ==  "Côte d’Ivoire"] <- "Ivory Coast"
  
  world_cluster <- world %>%
    left_join(cluster_map, by = c("region" = "country"))
  
  world_cluster$cluster <- as.character(world_cluster$cluster)
  ##############################################################################
  ## dont remove na, it will remove regions from the map, like greenland,
  
  label_df <- world_cluster %>%
    filter(region %in% cluster_map$country) %>%
    group_by(region) %>%
    summarise(
      long = mean(range(long)), # centroid longitude
      lat  = mean(range(lat)),  # centroid latitude
      cluster = first(cluster)
    )
  
  label_df$lat[label_df$region == "USA"]  <- 41.5
  label_df$long[label_df$region == "USA"] <- -100.00
  
  ##############################################################################
  
  wc_plot = ggplot(world_cluster,
                   aes(x = long, y = lat, group = group)) +
    geom_polygon(
      aes(fill = factor(cluster, levels = legend_levels)),
      color = "gray90",
      size = 0.01
    ) +
    geom_text(data = label_df, aes(x = long, y = lat, label = region),
              inherit.aes = FALSE, size = 1)+
    
    scale_fill_manual(
      values = cluster_colors,
      labels = legend_labels,
      name = "Clusters", drop   = FALSE
      
    ) +
    
    theme_void()+
    
    labs(
      title = paste0("World Clusters: ", input_year)
    )
  
  
  return(list(clustering_plot, wc_plot))
}
################# plot generation ##############################
library(future)
library(future.apply)
library(parallel)

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)

years <- 1990:2000

############### HETCOV PLOTS - with cluster nos =3,4,5 #########################

################# only gdp ###############################################
wdi_ind_list <- c(
  "NY.GDP.PCAP.CD"
)


out <- future_lapply(years, function(yr) {
  cluster_map_cl3 = het_cov_unga_final(input_year = yr, 
                                       K=3,
                                       method = c('het_yes_cov'), 
                                       wdi_indicators = c(
                                         "NY.GDP.PCAP.CD"
                                       ), 
                                       include_hdi = 'no')
  cluster_map_cl4 = het_cov_unga_final(input_year = yr, 
                                       K=4,
                                       method = c('het_yes_cov'), 
                                       wdi_indicators = c(
                                         "NY.GDP.PCAP.CD"
                                       ), 
                                       include_hdi = 'no')
  cluster_map_cl5 = het_cov_unga_final(input_year = yr, 
                                       K=5,
                                       method = c('het_yes_cov'), 
                                       wdi_indicators = c(
                                         "NY.GDP.PCAP.CD"
                                       ), 
                                       include_hdi = 'no')
  
  cl_plot3 = unga_plots(cluster_map_cl3, yr)[[1]]
  cl_plot4 = unga_plots(cluster_map_cl4, yr)[[1]]
  cl_plot5 = unga_plots(cluster_map_cl5, yr)[[1]]
  
  
  ### save at desired location ###
  path = "Plots/wc"
  
  ggsave(
    filename = paste0("clplot_hetcov_norm_cl=3_", yr, ".pdf"),
    plot     = cl_plot3, 
    path     = path,
    width    = 12,
    height   = 8,
    
  )
  
  ggsave(
    filename = paste0("clplot_hetcov_norm_cl=4_", yr, ".pdf"),
    plot     = cl_plot4, 
    path     = path,
    width    = 12,
    height   = 8,
    
  )
  
  ggsave(
    filename = paste0("clplot_hetcov_norm_cl=5_", yr, ".pdf"),
    plot     = cl_plot5, 
    path     = path,
    width    = 12,
    height   = 8,
    
  )
  
  
})

plan(sequential)


############### HETCOV PLOTS (cl_plot_cl=3,wc_plots)############################

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)


years <- 1990:2000

################# only gdp ###############################################
wdi_ind_list <- c(
  "NY.GDP.PCAP.CD"
)


out <- future_lapply(years, function(yr) {
  cluster_map = het_cov_unga_final(input_year = yr, 
                                   K=3,
                                   method = c('het_yes_cov'), 
                                   wdi_indicators = c(
                                     "NY.GDP.PCAP.CD"
                                   ), 
                                   include_hdi = 'no')
  
  cl_plot = unga_plots(cluster_map, yr)[[1]]
  
  wc_plot = unga_plots(cluster_map, yr)[[2]]
  
  ### save at desired location ###
  path = "Plots/wc"
  
  ggsave(
    filename = paste0("clplot_hetcov_norm_cl=3_", yr, ".pdf"),
    plot     = cl_plot, 
    path     = path,
    width    = 12,
    height   = 8,
    
  )
  
  ggsave(
    filename = paste0("wc_hetcov_norm_cl=3_", yr, ".pdf"),
    plot     = wc_plot,   
    path     = path,
    width    = 12,
    height   = 8,
    
  )
  
})



################### HET PLOTS ##################################################
out_het_cl_3 <- future_lapply(years, function(yr) {
  cluster_map = het_cov_unga_final(input_year = yr, 
                                   K=3,
                                   method = c('het_no_cov'), 
                                   wdi_indicators = c(
                                     "NY.GDP.PCAP.CD"), 
                                   include_hdi = 'no')
  
  cl_plot = unga_plots(cluster_map, yr)[[1]]
  
  wc_plot = unga_plots(cluster_map, yr)[[2]]
  
  ### save at desired location ###
  
  ggsave(
    filename = paste0("clplot_het_norm_cl=3_", yr, ".pdf"),
    plot     = cl_plot, 
    path     = "Plots/wc",
    width    = 12,
    height   = 8,
    
  )
  
  ggsave(
    filename = paste0("wc_het_norm_cl=3_", yr, ".pdf"),
    plot     = wc_plot,   
    path     = "Plots/wc",
    width    = 12,
    height   = 8,
    
  )
  
})

plan(sequential)

