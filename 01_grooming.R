library(cmdstanr)
pathtoresults <- "replication_pack/results/pan_grooming"

xdata <- read.csv("replication_pack/data_files/data_grooming.csv")
pdata <- read.csv("replication_pack/data_files/data_grooming_predictions.csv")

# prepare data for Stan
sdat <- list(n_obs = nrow(xdata),
             obseff = xdata$obs_hours,
             spec_index = (xdata$species == "bonobo") + 1,
             sexspec_indicator = xdata$sexspec_indicator,
             partners_obs = xdata$partners_groom,
             ntrials = xdata$partners_pool,
             index_indi = as.numeric(as.factor(xdata$individual)),
             n_indi = length(unique(xdata$individual)),
             groupsize = xdata$partners_pool + 1,
             n_groups = length(unique(xdata$groupyear)),
             index_group = as.numeric(as.factor(xdata$groupyear))
)
sdat$groupsize <- sdat$groupsize - mean(sdat$groupsize)
names(sdat$index_group) <- xdata$groupyear
names(sdat$sexspec_indicator) <- xdata$sexspec

lapply(sdat, head)


# prediction stuff
sdat$pred_size <- pdata$groupsize # group size
sdat$pred_sexspec_index <- pdata$pred_sexspec_index
names(sdat$pred_sexspec_index) <- pdata$pred_sexspec_index
sdat$pred_group_index <- pdata$pred_group_index
names(sdat$pred_group_index) <- pdata$group_year
sdat$pred_obseff <- 3 # target obseff for predictions
sdat$n_preds <- nrow(pdata)
lapply(sdat, head)






## fit model ----
m <- cmdstan_model("replication_pack/stan_files/pan_grooming.stan")
m$check_syntax(pedantic = TRUE)

# fit model
r <- m$sample(data = sdat, refresh = 0, parallel_chains = 4, seed = 1, init = 0.1)


rs <- r$summary(c("b", "indi_sd", "group_sd", "gs_slope"))
aux <- apply(r$draws(c("b", "indi_sd", "group_sd", "gs_slope"), format = "draws_matrix"), 2, quantile, probs = c(0.055, 0.945))
rs$ql89 <- aux[1, ]
rs$qu89 <- aux[2, ]
aux <- apply(r$draws(c("b", "indi_sd", "group_sd", "gs_slope"), format = "draws_matrix"), 2, quantile, probs = c(0.25, 0.75))
rs$ql50 <- aux[1, ]
rs$qu50 <- aux[2, ]


if (!dir.exists(pathtoresults)) dir.create(pathtoresults)
r$save_object(file.path(pathtoresults, "modelenv_pan_grooming.rds"))
saveRDS(rs, file = file.path(pathtoresults, "summarydataframe_pan_grooming.rds"))
saveRDS(sdat, file = file.path(pathtoresults, "standat_pan_grooming.rds"))
saveRDS(pdata, file = file.path(pathtoresults, "prediction_skeleton_pan_grooming.rds"))
saveRDS(xdata, file = file.path(pathtoresults, "rawdata_pan_grooming.rds"))

