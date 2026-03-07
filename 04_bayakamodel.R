# Bayaka cooperation model
pathtoresults <- "replication_pack/results/bayaka"
source("replication_pack/functions/prep_data.R")

cmdstanr::cmdstan_path()
library(cmdstanr)
m <- cmdstanr::cmdstan_model("replication_pack/stan_files/bayaka_model.stan")
m$check_syntax(pedantic = TRUE)

# read full data set
xdata <- read.csv("replication_pack/data_files/human_data.csv")
sdat <- prep_dat_hum(xdata = xdata)

# demo/debugging version with subset
set.seed(123)
xdata <- xdata[sample(nrow(xdata), nrow(xdata)/4), ]
sdat <- prep_dat_hum(xdata = xdata)
r <- m$sample(sdat, parallel_chains = 3, chains = 3, iter_warmup = 200, iter_sampling = 142, refresh = 0, show_exceptions = FALSE, seed = 42)
# takes 20 seconds on 5-year old laptop

# full set
# xdata <- read.csv("replication_pack/data_files/human_data.csv")
# sdat <- prep_dat_hum(xdata = xdata)
# r <- m$sample(data = sdat, parallel_chains = 6, chains = 12, refresh = 50,
#               iter_warmup = 1500, iter_sampling = 1000,
#               show_exceptions = FALSE,
#               max_treedepth = 13, adapt_delta = 0.85, seed = 1)


rs <- data.frame(r$summary()) # in demo model: warnings due to small number of samples
temp <- data.frame(r$summary(quantiles = ~ quantile(., probs = c(0.055, 0.945))))[, 2:3]
colnames(temp) <- c("ql89", "qu89")
rs <- cbind(rs, temp)

if (!dir.exists(pathtoresults)) dir.create(pathtoresults)
r$save_object(file.path(pathtoresults, "modelenv_bayaka_cooperation.rds"))
saveRDS(rs, file = file.path(pathtoresults, "summarydataframe_bayaka_cooperation.rds"))
saveRDS(sdat, file = file.path(pathtoresults, "standat_bayaka_cooperation.rds"))
