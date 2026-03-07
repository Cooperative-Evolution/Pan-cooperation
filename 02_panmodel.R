# Pan cooperation model
pathtoresults <- "replication_pack/results/pan_cooperation"
source("replication_pack/functions/prep_data.R")


cmdstanr::cmdstan_path()
panmod <- cmdstanr::cmdstan_model("replication_pack/stan_files/pan_cooperation.stan")



xdata <- read.csv("replication_pack/data_files/pandata.csv")

# toggle subset for illustration/debugging ----
# - demo1 5 or 6 individuals per group-year
# - demo2 7 or 8 individuals per group-year
xdata <- xdata[xdata$demo1, ]


# prep data set for model fit ----
xx <- prep_dat(xdata)
sdat <- xx$standat
ids <- xx$id_table
dyads <- xx$dyad_table
# peek
lapply(sdat, head)

# fit model ----
# full sample (takes 5 hours on a MacBook M1)
# careful: running on 8 cores in parallel (if you don't have 8: reduce parallel chains to 4 (this will double the time...))
# r <- panmod$sample(data = sdat, parallel_chains = 8, chains = 8, refresh = 10,
#                 iter_warmup = 1500, iter_sampling = 1000,
#                 max_treedepth = 13, adapt_delta = 0.99,
#                 seed = 47, show_exceptions = FALSE)
# debugging version
r <- panmod$sample(data = sdat, parallel_chains = 5, chains = 5, refresh = 50,
                   iter_warmup = 400, iter_sampling = 100, adapt_delta = 0.8,
                   seed = 754742, show_exceptions = FALSE, max_treedepth = 11)


rs <- r$summary(c("sigma_sexspec", "icpts_sexspec", "sigma_combispec", "icpts_combispec", "group_sd", "pop_sd", "year_sd", "b", "b_kin"))

if (!dir.exists(pathtoresults)) dir.create(pathtoresults)
r$save_object(file.path(pathtoresults, "modelenv_pan_cooperation.rds"))
saveRDS(rs, file = file.path(pathtoresults, "summarydataframe_pan_cooperation.rds"))
saveRDS(sdat, file = file.path(pathtoresults, "standat_pan_cooperation.rds"))


r$diagnostic_summary()

r$summary("icpts_sexspec")
