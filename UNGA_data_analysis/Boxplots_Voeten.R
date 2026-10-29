
library(dplyr); library(tidyr); library(ggplot2)
library(countrycode)   # for name harmonization
library(ggrepel)       # for non-overlapping labels
library(mclust)        # for ARI

### load het_cov_unga_final
source("Helpers.R")


######## BOXPLOTS ########################

compare_clusters_voeten <- function(input_year, K=3) {
  
  # ---- Load Voeten ideal points ----
  if (!file.exists("voeten_ip.tab")) {
    url <- "https://dataverse.harvard.edu/api/access/datafile/13642025"
    download.file(url, "voeten_ip.tab", mode = "wb")
  }
  ip <- read.delim("voeten_ip.tab", stringsAsFactors = FALSE)
  
  ip_year <- ip %>%
    filter(year == input_year) %>%
    transmute(country = Countryname,
              ideal_point = IdealPointFP) %>%
    mutate(iso3 = countrycode(country, "country.name", "iso3c"))
  
  # ---- Run hetcov clustering ----
  cl <- het_cov_unga_final(input_year = input_year,
                           K = K,
                           method = c('het_yes_cov'),
                           wdi_indicators = 'NY.GDP.PCAP.CD',
                           include_hdi = 'no')
  
  cl <- cl %>% mutate(iso3 = countrycode(country, "country.name", "iso3c"))
  
  # ---- Tag blocs based on anchor countries ----
  us_cl    <- cl$cluster[cl$country == "United States"]
  china_cl <- cl$cluster[cl$country == "China"]
  gulf_cl  <- cl$cluster[cl$country == "Qatar"]
  
  cl$bloc <- "Other"
  cl$bloc[cl$cluster == us_cl]    <- "USA bloc"
  cl$bloc[cl$cluster == china_cl] <- "China bloc"
  
  # Only assign Gulf bloc if K >= 4 AND Qatar is in its own cluster
  has_gulf <- K >= 4 &&
    length(gulf_cl) > 0 &&
    !(gulf_cl %in% c(us_cl, china_cl))
  
  if (has_gulf) {
    cl$bloc[cl$cluster == gulf_cl] <- "Gulf bloc"
  }
  
  # Set factor levels (drop Gulf bloc if not present)
  bloc_levels <- c("USA bloc", "China bloc")
  if (has_gulf) bloc_levels <- c(bloc_levels, "Gulf bloc")
  bloc_levels <- c(bloc_levels, "Other")
  cl$bloc <- factor(cl$bloc, levels = bloc_levels)
  
  # ---- Merge ----
  merged <- inner_join(cl, ip_year, by = "iso3") %>%
    mutate(country = country.x) %>%
    select(country, everything(), -country.x, -country.y)
  
  # ---- Palette (only keep the blocs that exist) ----
  pal_full <- c("USA bloc"   = "royalblue",
                "China bloc" = "#d7301f",
                "Gulf bloc"  = "#A65CA4",
                "Other"      = "#31a354")
  pal <- pal_full[bloc_levels]
  
  important <- c("United States", "China", "Russia", "Russian Federation","Portugal","Spain",
                 "India", "United Kingdom", "France", "Israel", "Cuba",
                 "North Korea", "Poland", "Hungary", "Germany", "Japan", "Iran",
                 "United Arab Emirates", "Qatar", "Bahrain", "Kuwait", "Oman", #gulf countries
                 "Argentina","Chile","Brazil","Bolivia","Paraguay","Venezuela","Mexico","Libya") # latin american countries
  
  # ---- Plot A: box + jitter ----
  boxplot <- ggplot(merged, aes(x = bloc, y = ideal_point, fill = bloc)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    geom_text_repel(data = filter(merged, country %in% important),
                    aes(label = country), size = 3) +
    scale_fill_manual(values = pal) +
    labs(title = paste("Ideal points by cluster —", input_year),
         x = NULL, y = "Ideal point")
  
  
  list(boxplot = boxplot)
  
}


###################PLOTS#############################################################
library(future)
library(future.apply)
library(parallel)

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)


years <- 1990:2000


future_lapply(years, function(yr) {
  
  res =compare_clusters_voeten(input_year = yr, K = 3)
  bp = res$boxplot
  
  
  ggsave(
    filename = paste0("bp_", yr, ".pdf"),
    plot     = bp, 
    path     = "Plots/boxplots_voeten",
    width    = 12,
    height   = 8,
    
  )
  
  
})

plan(sequential)
###############################################################################
