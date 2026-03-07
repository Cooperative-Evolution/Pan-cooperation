library(cmdstanr)
pathtoresults <- "replication_pack/results/pan_cooperation/"
source("replication_pack/functions/plot_post_cors.R")

# reread results from disk ----
r <- readRDS(file.path(pathtoresults, "modelenv_pan_cooperation.rds"))
sdat <- readRDS(file.path(pathtoresults, "standat_pan_cooperation.rds"))
rs <- readRDS(file.path(pathtoresults, "summarydataframe_pan_cooperation.rds"))


# ______________________________________________________________________________
# figure 2 ---- correlations of dyad-specific propensities
# ______________________________________________________________________________

plot_post_cor_4axes(model_env = r, standat = sdat, by_individual = FALSE)
