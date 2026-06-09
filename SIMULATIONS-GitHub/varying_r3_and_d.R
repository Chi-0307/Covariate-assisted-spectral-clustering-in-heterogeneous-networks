
######===== Descritption of the code =============#######
# 1.Vary r3 and dcov to plot the misclustering error rate for all the
# baselines and hetcov.
# 2.Joint effect of r3 and dcov is also plotted.

################################################################################
library(ggplot2);library(dplyr);library(tidyr)
library(future)
library(future.apply)
###############################################################################
source("Helper functions.R")
################################################################################
#################### vary r3 - faster code #####################################

## parallel plan (Windows-safe)
future::plan(future::multisession, workers = future::availableCores() - 1)

## parameters
r3_vals <- seq(0, 0.7, 0.05)
nsim <- 50
n <- 250

## storage
n_hetcov1 <- matrix(NA, length(r3_vals), nsim)
n_hetcov2 <- matrix(NA, length(r3_vals), nsim)

n_het1 <- matrix(NA, length(r3_vals), nsim)
n_het2 <- matrix(NA, length(r3_vals), nsim)

n_homcov1 <- matrix(NA, length(r3_vals), nsim)
n_homcov2 <- matrix(NA, length(r3_vals), nsim)

n_hom1 <- matrix(NA, length(r3_vals), nsim)
n_hom2 <- matrix(NA, length(r3_vals), nsim)

misclassification_rate_hetcov_r3_1 <- numeric(length(r3_vals))
misclassification_rate_hetcov_r3_2 <- numeric(length(r3_vals))

misclassification_rate_het_r3_1 <- numeric(length(r3_vals))
misclassification_rate_het_r3_2 <- numeric(length(r3_vals))

misclassification_rate_homcov_r3_1 <- numeric(length(r3_vals))
misclassification_rate_homcov_r3_2 <- numeric(length(r3_vals))

misclassification_rate_hom_r3_1 <- numeric(length(r3_vals))
misclassification_rate_hom_r3_2 <- numeric(length(r3_vals))

## reproducibility
set.seed(1234)

## main loop
for (i in seq_along(r3_vals)) {
  
  r3 <- r3_vals[i]
  
  res <- future.apply::future_lapply(
    1:nsim,
    function(sim) {
      #set.seed(123)
      
      list(                               
        hetcov1 = het_cov_sbm(
          d = 0.5, r = 0.01, r3 = r3,
          n = n, K = 3, T = 2, R = 3
        ),
        
        hetcov2 = het_cov_sbm(
          d = 0.5, r = 0.01, r3 = r3,
          n = n, K = 3, T = 2, R = 3
        ),
        
        het1 = heterogeneous_sbm(
          r = 0.01, r3 = r3,
          n = n, K = 3, T = 2
        ),
        
        het2 = heterogeneous_sbm(
          r = 0.01, r3 = r3,
          n = n, K = 3, T = 2
        ),
        
        homcov1 = homogeneous_cov_sbm(
          d = 0.5, r = 0.01,
          n = n, K = 3, T = 2, R = 3
        ),
        
        homcov2 = homogeneous_cov_sbm(
          d = 0.5, r = 0.01,
          n = n, K = 3, T = 2, R = 3
        ),
        
        hom1 = homogeneous_sbm(
          r = 0.01,
          n = n, K = 3
        ),
        
        hom2 = homogeneous_sbm(
          r = 0.01,
          n = n, K = 3
        )
      )
    },
    future.seed = TRUE
  )
  
  ## extract
  n_hetcov1[i, ] <- sapply(res, function(x) x$hetcov1[[1]])
  n_hetcov2[i, ] <- sapply(res, function(x) x$hetcov2[[2]])
  
  n_het1[i, ] <- sapply(res, function(x) x$het1[[1]])
  n_het2[i, ] <- sapply(res, function(x) x$het2[[2]])
  
  n_homcov1[i, ] <- sapply(res, function(x) x$homcov1[[1]])
  n_homcov2[i, ] <- sapply(res, function(x) x$homcov2[[2]])
  
  n_hom1[i, ] <- sapply(res, function(x) x$hom1[[1]])
  n_hom2[i, ] <- sapply(res, function(x) x$hom2[[2]])
  
  ## misclassification rates
  misclassification_rate_hetcov_r3_1[i] <- mean(n_hetcov1[i, ]) / (n / 2)
  misclassification_rate_hetcov_r3_2[i] <- mean(n_hetcov2[i, ]) / (n / 2)
  
  misclassification_rate_het_r3_1[i] <- mean(n_het1[i, ]) / (n / 2)
  misclassification_rate_het_r3_2[i] <- mean(n_het2[i, ]) / (n / 2)
  
  misclassification_rate_homcov_r3_1[i] <- mean(n_homcov1[i, ]) / (n / 2)
  misclassification_rate_homcov_r3_2[i] <- mean(n_homcov2[i, ]) / (n / 2)
  
  misclassification_rate_hom_r3_1[i] <- mean(n_hom1[i, ]) / (n / 2)
  misclassification_rate_hom_r3_2[i] <- mean(n_hom2[i, ]) / (n / 2)
}

## shutdown
future::plan(sequential)


############## PLOT r3 #########################################################
library(ggplot2)

# Combine both panels into one dataframe
df <- data.frame(
  r3 = rep(r3_vals, 8),
  value = c(
    misclassification_rate_hetcov_r3_1,
    misclassification_rate_het_r3_1,
    misclassification_rate_homcov_r3_1,
    misclassification_rate_hom_r3_1,
    misclassification_rate_hetcov_r3_2,
    misclassification_rate_het_r3_2,
    misclassification_rate_homcov_r3_2,
    misclassification_rate_hom_r3_2
  ),
  method = factor(rep(
    c("het-cov", "het", "hom-cov", "hom"),
    each = length(r3_vals),
    times = 2
  )),
  panel = factor(rep(c("Type 1", "Type 2"), each = 4 * length(r3_vals))) 
)

# Plot WITHOUT CONFIDENCE BANDS
ggplot(df, aes(x = r3, y = value, color = method)) +
  geom_line(size = 0.8) + 
  facet_wrap(~panel, nrow = 1) +   # side-by-side
  scale_color_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    breaks = c("het-cov", "het", "hom-cov", "hom"),  # controls order in legend
    name = "Method"
  ) +
  labs(
    x = expression(r[3]),
    y = "Average Misclassification Error Rate"
  ) +
  ylim(0, 0.8) +
  theme_minimal() 


############ mean +- sd plots ##################################

get_stats <- function(mat, n) {
  data.frame(
    mean = rowMeans(mat) / n,
    sd   = apply(mat, 1, sd) / n
  )
}
# no of type1 nodes
n1 = n/2
# no of type2 nodes
n2 = n1

# Build combined dataframe
df_all <- data.frame(
  r3 = rep(r3_vals, 8),
  
  mean = c(
    get_stats(n_hetcov1, n1)$mean,
    get_stats(n_het1, n1)$mean,
    get_stats(n_homcov1, n1)$mean,
    get_stats(n_hom1, n1)$mean,
    get_stats(n_hetcov2, n2)$mean,
    get_stats(n_het2, n2)$mean,
    get_stats(n_homcov2, n2)$mean,
    get_stats(n_hom2, n2)$mean
  ),
  
  sd = c(
    get_stats(n_hetcov1, n1)$sd,
    get_stats(n_het1, n1)$sd,
    get_stats(n_homcov1, n1)$sd,
    get_stats(n_hom1, n1)$sd,
    get_stats(n_hetcov2, n2)$sd,
    get_stats(n_het2, n2)$sd,
    get_stats(n_homcov2, n2)$sd,
    get_stats(n_hom2, n2)$sd
  ),
  
  method = factor(rep(
    c("het-cov", "het", "hom-cov", "hom"),
    each = length(r3_vals),
    levels = c("het-cov", "het", "hom-cov", "hom"),
    times = 2
  )),
  
  panel = factor(rep(
    c("Type 1", "Type 2"),
    each = 4 * length(r3_vals)
  ))
)

df_all$method <- factor(
  df_all$method,
  levels = c("het-cov", "het", "hom-cov", "hom")
)

# Plot
ggplot(df_all, aes(x = r3, y = mean, color = method, fill = method)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              alpha = 0.2, color = NA) +
  facet_wrap(~panel, nrow = 1) +
  scale_color_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    breaks = c("het-cov", "het", "hom-cov", "hom"),
    labels = c("Het-Cov", "Het", "Hom-Cov", "Hom"),
    name = "Method"
  )+
  scale_fill_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    guide = "none"
  ) +
  labs(
    x = expression(r[3]),
    y = "Average Misclassification Error Rate"
  ) + 
  ylim(0, 1) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 8)
  )



######################### vary d - FASTER CODE ##############################################################


library(future)
library(future.apply)

## parallel plan (Windows-safe)
future::plan(future::multisession, workers = future::availableCores() - 1)

## parameters
d_vals <- seq(0, 0.94, 0.05)
nsim <- 50
n <- 250
  
## storage
n_hetcov1 <- matrix(NA, length(d_vals), nsim)
n_hetcov2 <- matrix(NA, length(d_vals), nsim)

n_het1 <- matrix(NA, length(d_vals), nsim)
n_het2 <- matrix(NA, length(d_vals), nsim)

n_homcov1 <- matrix(NA, length(d_vals), nsim)
n_homcov2 <- matrix(NA, length(d_vals), nsim)

n_hom1 <- matrix(NA, length(d_vals), nsim)
n_hom2 <- matrix(NA, length(d_vals), nsim)

misclassification_rate_hetcov_d_1 <- numeric(length(d_vals))
misclassification_rate_hetcov_d_2 <- numeric(length(d_vals))

misclassification_rate_het_d_1 <- numeric(length(d_vals))
misclassification_rate_het_d_2 <- numeric(length(d_vals))

misclassification_rate_homcov_d_1 <- numeric(length(d_vals))
misclassification_rate_homcov_d_2 <- numeric(length(d_vals))

misclassification_rate_hom_d_1 <- numeric(length(d_vals))
misclassification_rate_hom_d_2 <- numeric(length(d_vals))

## reproducibility
set.seed(1234)

## main loop
for (i in seq_along(d_vals)) {
  
  d <- d_vals[i]
  
  res <- future.apply::future_lapply(
    1:nsim,
    function(sim) {
      #set.seed(123)
      
      list(
        hetcov1 = het_cov_sbm(
          d = d, r = 0.01, r3 = 0.2,
          n = n, K = 3, T = 2, R = 3
        ),
        
        hetcov2 = het_cov_sbm(
          d = d, r = 0.01, r3 = 0.2,
          n = n, K = 3, T = 2, R = 3
        ),
        
        het1 = heterogeneous_sbm(
          r = 0.01, r3 = 0.2,
          n = n, K = 3, T = 2
        ),
        
        het2 = heterogeneous_sbm(
          r = 0.01, r3 = 0.2,
          n = n, K = 3, T = 2
        ),
        
        homcov1 = homogeneous_cov_sbm(
          d = d, r = 0.01,
          n = n, K = 3, T = 2, R = 3
        ),
        
        homcov2 = homogeneous_cov_sbm(
          d = d, r = 0.01,
          n = n, K = 3, T = 2, R = 3
        ),
        
        hom1 = homogeneous_sbm(
          r = 0.01,
          n = n, K = 3
        ),
        
        hom2 = homogeneous_sbm(
          r = 0.01,
          n = n, K = 3
        )
      )
    },
    future.seed = TRUE
  )
  
  ## extract
  n_hetcov1[i, ] <- sapply(res, function(x) x$hetcov1[[1]])
  n_hetcov2[i, ] <- sapply(res, function(x) x$hetcov2[[2]])
  
  n_het1[i, ] <- sapply(res, function(x) x$het1[[1]])
  n_het2[i, ] <- sapply(res, function(x) x$het2[[2]])
  
  n_homcov1[i, ] <- sapply(res, function(x) x$homcov1[[1]])
  n_homcov2[i, ] <- sapply(res, function(x) x$homcov2[[2]])
  
  n_hom1[i, ] <- sapply(res, function(x) x$hom1[[1]])
  n_hom2[i, ] <- sapply(res, function(x) x$hom2[[2]])
  
  ## misclassification rates
  misclassification_rate_hetcov_d_1[i] <- mean(n_hetcov1[i, ]) / (n / 2)
  misclassification_rate_hetcov_d_2[i] <- mean(n_hetcov2[i, ]) / (n / 2)
  
  misclassification_rate_het_d_1[i] <- mean(n_het1[i, ]) / (n / 2)
  misclassification_rate_het_d_2[i] <- mean(n_het2[i, ]) / (n / 2)
  
  misclassification_rate_homcov_d_1[i] <- mean(n_homcov1[i, ]) / (n / 2)
  misclassification_rate_homcov_d_2[i] <- mean(n_homcov2[i, ]) / (n / 2)
  
  misclassification_rate_hom_d_1[i] <- mean(n_hom1[i, ]) / (n / 2)
  misclassification_rate_hom_d_2[i] <- mean(n_hom2[i, ]) / (n / 2)
}

## shutdown
future::plan(sequential)

################ PLOT d - WITHOUT CI ##########################################################
df <- data.frame(
  d = rep(d_vals, 8),
  value = c(
    misclassification_rate_hetcov_d_1,
    misclassification_rate_het_d_1,
    misclassification_rate_homcov_d_1,
    misclassification_rate_hom_d_1,
    misclassification_rate_hetcov_d_2,
    misclassification_rate_het_d_2,
    misclassification_rate_homcov_d_2,
    misclassification_rate_hom_d_2
  ),
  method = factor(rep(
    c("het-cov", "het", "hom-cov", "hom"),
    each = length(d_vals),
    times = 2
  )),
  panel = factor(rep(c("Type 1", "Type 2"), each = 4 * length(d_vals))) 
)

# Plot with facets
ggplot(df, aes(x = d, y = value, color = method)) +
  geom_line(size = 0.5) + 
  facet_wrap(~panel, nrow = 1) +   # side-by-side
  scale_color_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    breaks = c("het-cov", "het", "hom-cov", "hom"),  # controls order in legend
    name = "Method"
  ) +
  labs(
    x = expression(d[cov]),
    y = "Average Misclassification Error Rate"
  ) +
  ylim(0, 1) +
  theme_minimal() 

################################################################################
################ PLOT d (with CI)################################################

get_stats <- function(mat, n) {
  data.frame(
    mean = rowMeans(mat) / (n/2),
    sd   = apply(mat, 1, sd) / (n/2)
  )
}

df1 <- rbind(
  cbind(d = d_vals, get_stats(n_hetcov1, n), method = "het-cov", panel = "Type 1"),
  cbind(d = d_vals, get_stats(n_het1, n), method = "het", panel = "Type 1"),
  cbind(d = d_vals, get_stats(n_homcov1, n), method = "hom-cov", panel = "Type 1"),
  cbind(d = d_vals, get_stats(n_hom1, n), method = "hom", panel = "Type 1"),
  
  cbind(d = d_vals, get_stats(n_hetcov2, n), method = "het-cov", panel = "Type 2"),
  cbind(d = d_vals, get_stats(n_het2, n), method = "het", panel = "Type 2"),
  cbind(d = d_vals, get_stats(n_homcov2, n), method = "hom-cov", panel = "Type 2"),
  cbind(d = d_vals, get_stats(n_hom2, n), method = "hom", panel = "Type 2")
)

ggplot(df1, aes(x = d, y = mean, color = method, fill = method)) +
  geom_line(linewidth = 0.8) +
  
  geom_ribbon(
    aes(ymin = mean - sd, ymax = mean + sd),
    alpha = 0.2,
    color = NA
  ) +
  
  facet_wrap(~panel, nrow = 1) +
  
  scale_color_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    breaks = c("het-cov", "het", "hom-cov", "hom"),
    labels = c("Het-Cov", "Het", "Hom-Cov", "Hom"),
    name = "Method"
  ) +
  
  scale_fill_manual(
    values = c(
      "het-cov" = "firebrick1",
      "het" = "royalblue",
      "hom-cov" = "forestgreen",
      "hom" = "black"
    ),
    guide = "none"
  ) +
  
  labs(
    x = expression(d[cov]),
    y = "Average Misclassification Error Rate"
  ) +
  
  ylim(0, 1) +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 8)
  )

plan(sequential)

################################################################################
############### vary both (d,r3) - faster code #################################

## PARAMETERS
r3_vals = seq(0, 0.55, 0.05)
d_vals  = seq(0, 0.55, 0.05)
nsim <- 50

##########################################
## CREATE DIAGONAL (d, r3) PAIRS -> 12 points
tuples <- Map(function(d, r3) c(d, r3), d_vals, r3_vals)

## USE ALL POINTS
ticks <- tuples   # length = 12

##########################################
## ----------- PLOTTING ----------- ##

x_index <- 1:length(ticks)

## labels as (r3, d)
tick_labels <- sapply(ticks, function(x) {
  paste0("(", round(x[2], 2), ",", round(x[1], 2), ")")
})


## parallel plan
#future::plan(future::multisession, workers = 2)
future::plan(future::multisession, workers = future::availableCores() - 1)


##########################################
# type-1
n_hetcov_n_250_varyboth_t1 <- matrix(NA, length(ticks), nsim)
n_hetcov_n_500_varyboth_t1 <- matrix(NA, length(ticks), nsim)

## misclassification rate
misclassification_rate_hetcov_n_250_varyboth_t1 <- numeric(length(ticks))
misclassification_rate_hetcov_n_500_varyboth_t1 <- numeric(length(ticks))

# type-2
n_hetcov_n_250_varyboth_t2 <- matrix(NA, length(ticks), nsim)
n_hetcov_n_500_varyboth_t2 <- matrix(NA, length(ticks), nsim)

## misclassification rate
misclassification_rate_hetcov_n_250_varyboth_t2 <- numeric(length(ticks))
misclassification_rate_hetcov_n_500_varyboth_t2 <- numeric(length(ticks))

## reproducibility
set.seed(1234)

## main loop
for (i in seq_along(ticks)) {
  d <- ticks[[i]][1]
  r3 <- ticks[[i]][2]
  res <- future.apply::future_lapply(
    1:nsim,
    function(sim) {
      #set.seed(123)
      list(  
        #### n = 250 ###
        hetcov_n_250 = het_cov_sbm(
          d = d, r = 0.01, r3 = r3,
          n = 250, K = 3, T = 2, R = 3
        ),
        
        #### n = 500 ###
        hetcov_n_500 = het_cov_sbm(
          d = d, r = 0.01, r3 = r3,
          n = 500, K = 3, T = 2, R = 3
        )
        
      )
    },
    future.seed = TRUE
  )
  
  ## extract
  ### n = 250 #####
  n = 250
  n_hetcov_n_250_varyboth_t1[i, ] <- sapply(res, function(x) x$hetcov_n_250[[1]])
  n_hetcov_n_250_varyboth_t2[i, ] <- sapply(res, function(x) x$hetcov_n_250[[2]])
  # miscl
  misclassification_rate_hetcov_n_250_varyboth_t1[i] <- mean(n_hetcov_n_250_varyboth_t1[i, ]) / (n / 2)
  misclassification_rate_hetcov_n_250_varyboth_t2[i] <- mean(n_hetcov_n_250_varyboth_t2[i, ]) / (n / 2)
  
  ### n = 500 ####
  n = 500
  n_hetcov_n_500_varyboth_t1[i, ] <- sapply(res, function(x) x$hetcov_n_500[[1]])
  n_hetcov_n_500_varyboth_t2[i, ] <- sapply(res, function(x) x$hetcov_n_500[[2]])
  # miscl
  misclassification_rate_hetcov_n_500_varyboth_t1[i] <- mean(n_hetcov_n_500_varyboth_t1[i, ]) / (n / 2)
  misclassification_rate_hetcov_n_500_varyboth_t2[i] <- mean(n_hetcov_n_500_varyboth_t2[i, ]) / (n / 2)
  
}

## shutdown
future::plan(sequential)



######################################################################
############## plot: vary both (d,r3) ################################
## PARAMETERS
r3_vals = seq(0, 0.55, 0.05)
d_vals  = seq(0, 0.55, 0.05)

##########################################
## CREATE DIAGONAL (d, r3) PAIRS - 12 points
tuples <- Map(function(d, r3) c(d, r3), d_vals, r3_vals)

## USE ALL POINTS
ticks <- tuples   # length = 12

##########################################
## ---- PLOTTING ----------- ##

x_index <- 1:length(ticks)

## labels as (r3, d)
tick_labels <- sapply(ticks, function(x) {
  paste0("(", round(x[2], 2), ",", round(x[1], 2), ")")
})

##########################################
## SIDE-BY-SIDE PLOTS
par(mfrow = c(1, 2), mar = c(7, 4, 3, 1), bg = "grey98")

##########################################
## TYPE-1
plot(
  x_index,
  misclassification_rate_hetcov_n_250_varyboth_t1,
  type = "l",
  col = "black",
  lwd = 1,
  ylim = c(0, 1),
  xaxt = "n",
  xlab = expression((r[3] * "," ~ d[cov])),
  ylab = "Average Misclustering Rate",
  cex.lab = 0.8,
  main = "Type-1",
  font.main = 1,     # NOT bold
  cex.main = 0.9     # slightly smaller
)

lines(
  x_index,
  misclassification_rate_hetcov_n_500_varyboth_t1,
  col = "red",
  lwd = 1,
)



text(
  x = x_index,
  y = par("usr")[3] - 0.05,   # place below axis
  labels = tick_labels,
  srt = 45,                   # 45° rotation
  adj = 1,
  xpd = TRUE,
  cex = 0.7
)

grid(col = "grey80", lwd = 1)

legend(
  "topright",
  legend = c(expression(n == 250), expression(n == 500)),
  col = c("black", "red"),
  lwd = 2,
  bty = "n",     # removes box
  cex = 0.8,
  y.intersp = 0.5
  
)

##########################################
## TYPE-2
plot(
  x_index,
  misclassification_rate_hetcov_n_250_varyboth_t2,
  type = "l",
  col = "black",
  lwd = 1,
  ylim = c(0, 1),
  xaxt = "n",
  xlab = expression((r[3] * "," ~ d[cov])),
  ylab = "Average Misclustering Rate",
  cex.lab = 0.8,    ## controls font of ylab
  main = "Type-2",
  font.main = 1,     # NOT bold
  cex.main = 0.9     # slightly smaller
)

lines(
  x_index,
  misclassification_rate_hetcov_n_500_varyboth_t2,
  col = "red",
  lwd = 1
)

# axis(
#   1,
#   at = x_index,
#   labels = tick_labels,
#   las = 2,
#   cex.axis = 0.6
# )

text(
  x = x_index,
  y = par("usr")[3] - 0.05,   # place below axis
  labels = tick_labels,
  srt = 45,                   # 45° rotation
  adj = 1,
  xpd = TRUE,
  cex = 0.7
)


grid(col = "grey80", lwd = 1)

legend(
  "topright",
  legend = c(expression(n == 250), expression(n == 500)),
  col = c("black", "red"),
  lwd = 2,
  bty = "n",     # removes box
  cex = 0.8,
  y.intersp = 0.5
  
)

##########################################
## RESET
par(mfrow = c(1,1))

plan(sequential)
