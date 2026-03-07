# Pan association model
pathtoresults <- "replication_pack/results/pan_associations"

# install bamoso package
# remotes::install_github("gobbios/bamoso", dependencies = TRUE, build_vignettes = FALSE)
library(cmdstanr)
library(bamoso)
source("replication_pack/functions/prep_data.R")

# data prep ----
xdata <- read.csv("replication_pack/data_files/pandata.csv")
xdata$groupyear <- paste(xdata$group, xdata$year)
xdata$ind1 <- xdata$ind1_ano
xdata$ind2 <- xdata$ind2_ano


iter_warmup <- 5000
chains <- 6
iter_sampling <- 5000

# settings for demo/debugging
xdata <- xdata[xdata$demo1, ]
iter_warmup <- 500
chains <- 3
iter_sampling <- 100




## loop through groups and extract relevant draws
resmat_intercept <- matrix(nrow = chains * iter_sampling, ncol = length(unique(xdata$groupyear)))
colnames(resmat_intercept) <- unique(xdata$groupyear)
resmat_indi <- resmat_intercept
resmat_dyad <- resmat_intercept
summarylist <- list()

prior_mats <- matrix(ncol = 2, nrow = 0)

i=sort(unique(xdata$groupyear))[1]
for (i in unique(xdata$groupyear)) {
  cat("running", i, "\n")
  d <- make_assoc_standat(xdata, groupyear = i)
  prior_mats <- rbind(prior_mats, d$prior_matrix)
  zz <- capture.output( # in a silent way
    m <- suppressMessages(sociality_model(d,
                                          parallel_chains = 6,
                                          chains = chains,
                                          refresh = 0,
                                          adapt_delta = 0.99,
                                          max_treedepth = 13,
                                          iter_warmup = iter_warmup,
                                          iter_sampling = iter_sampling,
                                          seed = 12))
  )

  temp <- data.frame(m$mod_res$summary(c("beh_intercepts", "indi_soc_sd", "dyad_soc_sd")))
  temp <- data.frame(groupyear = i, temp)
  summarylist[[length(summarylist) + 1]] <- temp
  resmat_intercept[, i] <- c(m$mod_res$draws(c("beh_intercepts"), format = "draws_matrix"))
  resmat_indi[, i] <- c(m$mod_res$draws(c("indi_soc_sd"), format = "draws_matrix"))
  resmat_dyad[, i] <- c(m$mod_res$draws(c("dyad_soc_sd"), format = "draws_matrix"))
  cat("---------------------\n\n")
}


saveRDS(resmat_intercept, file = file.path(pathtoresults, "postdraws_baseline.rds"))
saveRDS(resmat_indi, file = file.path(pathtoresults, "postdraws_indivar.rds"))
saveRDS(resmat_dyad, file = file.path(pathtoresults, "postdraws_dyadvar.rds"))

