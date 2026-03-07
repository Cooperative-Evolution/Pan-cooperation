prep_dat <- function(xdata) {
  if (interactive()) {
    message("convert table format in to list for Stan")
    message("also creates an individual-level and a dyad-level table for post-processing")
  }

  xdata$ind1 <- xdata$ind1_ano
  xdata$ind2 <- xdata$ind2_ano

  # prep the data for Stan ----
  ids <- data.frame(id = unique(c(xdata$ind1, xdata$ind2)))
  ids$num <- seq_len(nrow(ids))
  ids$spec <- NA
  ids$group <- NA
  rownames(ids) <- ids$id
  ids$sex <- NA
  i=ids$id[1]
  for (i in ids$id) {
    aux <- xdata[xdata$ind1 == i, c("ind1", "sex1", "species", "group")]
    if (nrow(aux) == 0) aux <- xdata[xdata$ind2 == i, c("ind2", "sex2", "species", "group")]
    if (aux[1, 2] == "M") ids[i, "sex"] <- "M"
    if (aux[1, 2] == "F") ids[i, "sex"] <- "F"
    ids[i, "spec"] <- aux[1, "species"]
    ids[i, "group"] <- aux[1, "group"]
  }

  ids$specsex <- NA
  ids$specsex[ids$spec == "bonobo" & ids$sex == "F"] <- 1
  ids$specsex[ids$spec == "bonobo" & ids$sex == "M"] <- 2
  ids$specsex[ids$spec == "chimp" & ids$sex == "F"] <- 3
  ids$specsex[ids$spec == "chimp" & ids$sex == "M"] <- 4


  xdata$i1 <- t(ids[, "num", drop = FALSE])[1, xdata$ind1]
  xdata$i2 <- t(ids[, "num", drop = FALSE])[1, xdata$ind2]

  xdata$dyad <- apply(xdata[, c("i1", "i2")], 1, function(x)paste(sort(x), collapse = "_"))
  xdata$dyad_char <- apply(xdata[, c("ind1", "ind2")], 1, function(x)paste(sort(x), collapse = "_"))

  dyads <- data.frame(dyad = unique(xdata$dyad), dyad_char = unique(xdata$dyad_char))
  dyads$num <- seq_len(nrow(dyads))
  rownames(dyads) <- dyads$dyad
  head(dyads)
  dyads$combi <- NA
  dyads$spec <- NA
  dyads$is_kin <- NA
  dyads$group <- NA
  d=dyads$dyad[2]
  for (d in dyads$dyad) {
    auxids <- as.numeric(unlist(strsplit(d, "_")))
    auxsex <- sort(as.character(ids$sex[ids$num %in% auxids]))
    auxspec <- ids$spec[ids$num %in% auxids]
    auxgroup <- ids$group[ids$num %in% auxids]
    if (length(unique(auxspec)) != 1) stop("sexcombination error")
    if (length(unique(auxgroup)) != 1) stop("group error")
    dyads$spec[dyads$dyad == d] <- auxspec[1]
    if (identical(auxsex, c("F", "F"))) dyads$combi[dyads$dyad == d] <- "ff"
    if (identical(auxsex, c("M", "M"))) dyads$combi[dyads$dyad == d] <- "mm"
    if (identical(auxsex, c("F", "M"))) dyads$combi[dyads$dyad == d] <- "mix"
    # kin?
    auxkin <- which((xdata$i1 == auxids[1] & xdata$i2 == auxids[2]) | (xdata$i1 == auxids[2] & xdata$i2 == auxids[1]))
    auxkin <- xdata$kin[auxkin]
    if (length(unique(auxkin)) != 1) stop("kin error")
    dyads[d, "is_kin"] <- auxkin[1]
    dyads[d, "group"] <- auxgroup[1]
  }
  table(dyads$combi, useNA = "a")
  table(dyads$spec, useNA = "a")
  table(dyads$is_kin, useNA = "a")
  if (any(is.na(dyads$combi))) stop("")
  if (any(is.na(dyads$spec))) stop("")
  if (any(is.na(dyads$is_kin))) stop("")

  dyads$combi_num <- 0
  dyads$combi_num[dyads$spec == "bonobo" & dyads$combi == "ff"] <- 1
  dyads$combi_num[dyads$spec == "bonobo" & dyads$combi == "mix"] <- 2
  dyads$combi_num[dyads$spec == "bonobo" & dyads$combi == "mm"] <- 3
  dyads$combi_num[dyads$spec == "chimp" & dyads$combi == "ff"] <- 4
  dyads$combi_num[dyads$spec == "chimp" & dyads$combi == "mix"] <- 5
  dyads$combi_num[dyads$spec == "chimp" & dyads$combi == "mm"] <- 6

  table(dyads$combi_num)


  xdata$dyad_num <- NA
  # xdata$combi <- NA
  for (i in seq_len(nrow(xdata))) xdata$dyad_num[i] <- dyads$num[dyads$dyad == xdata$dyad[i]]
  # for (i in seq_len(nrow(xdata))) xdata$combi[i] <- dyads$combi[dyads$dyad == xdata$dyad[i]]
  # table(xdata$sex1, xdata$combi)
  # table(xdata$sex2, xdata$combi)
  #
  all(rowSums(table(xdata$dyad_num, xdata$dyad) != 0) == 1)


  xdata$combi <- paste(xdata$species, xdata$species)
  xdata$combi_num <- as.integer(as.factor(xdata$combi))

  head(xdata)

  # interactions ----
  ## grooming ----
  # replace zeros in grooming with tiny non-zero values
  # do it on group-year level
  orirng <- .Random.seed
  set.seed(1)
  groups <- unique(xdata$group)
  xdata$groom_prop <- xdata$grooming
  xdata$obs.tot <- xdata$obs_hours * 3600
  g=groups[1]
  for (g in groups) {
    years <- unique(xdata$year[xdata$group == g])
    y=years[1]
    for (y in years) {
      gy_index <- which(xdata$group == g & xdata$year == y)
      aux <- xdata$groom_prop[gy_index] / xdata$obs.tot[gy_index]
      newmedian <- median(aux[aux > 0]) / 1000
      if (any(aux == 0)) {
        runif(sum(aux == 0), newmedian/10, newmedian/5) * xdata$obs.tot[gy_index][aux == 0]
        xdata$groom_prop[gy_index][aux == 0] <- runif(sum(aux == 0), newmedian/10, newmedian/5) * xdata$obs.tot[gy_index][aux == 0]
        # sanity check 2
        if (!max(xdata$groom_prop[gy_index][aux == 0]) < min(xdata$groom_prop[gy_index][aux > 0])) stop("failure at sanity check 2", call. = FALSE)
      }

    }
  }
  .Random.seed <- orirng

  xdata$groom_prop <- xdata$groom_prop/(xdata$obs_hours * 3600)

  standat <- list(
    n_obs = nrow(xdata),
    n_ids = nrow(ids),
    n_dyads = nrow(dyads),
    n_behs = 4,
    n_cors = 6,
    n_grps = length(unique(xdata$group)),
    n_years = length(unique(xdata$groupyear)),
    n_pops = length(unique(xdata$population)),
    id1_for_obs = xdata$i1,
    id2_for_obs = xdata$i2,
    dyad_for_obs = xdata$dyad_num,
    group_for_obs = xdata$group_num,
    year_for_obs = as.integer(as.factor(xdata$groupyear)),
    pop_for_obs = as.integer(as.factor(xdata$population)),
    specsex_for_id = ids$specsex,
    combi_for_dyad = dyads$combi_num,
    is_kin = dyads$is_kin,
    obseff = xdata$obs_hours/100, # PER 100 HOURS
    agg = xdata$aggression,
    coal = xdata$coalition,
    foodshare = xdata$sharing,
    obseff_adlib = xdata$association2/100, # PER 100 scans
    grooming = xdata$groom_prop,
    proxim = xdata$association2/xdata$association.total
  )

  # standat$n_cors <- standat$n_behs * (standat$n_behs - 1) / 2
  # if ("coal" %in% names(standat) & "foodshare" %in% names(standat)) {
  #   standat$n_behs <- 5
  #   standat$n_cors <- standat$n_behs * (standat$n_behs - 1) / 2
  # }

  lapply(standat, head)

  standat$sexspec1_ids = which(ids$specsex == 1)
  standat$n_sexspec1_ids = length(standat$sexspec1_ids)
  standat$sexspec2_ids = which(ids$specsex == 2)
  standat$n_sexspec2_ids = length(standat$sexspec2_ids)
  standat$sexspec3_ids = which(ids$specsex == 3)
  standat$n_sexspec3_ids = length(standat$sexspec3_ids)
  standat$sexspec4_ids = which(ids$specsex == 4)
  standat$n_sexspec4_ids = length(standat$sexspec4_ids)

  standat$combispec1_dyads = which(dyads$combi_num == 1)
  standat$n_combispec1_dyads = length(standat$combispec1_dyads)
  standat$combispec2_dyads = which(dyads$combi_num == 2)
  standat$n_combispec2_dyads = length(standat$combispec2_dyads)
  standat$combispec3_dyads = which(dyads$combi_num == 3)
  standat$n_combispec3_dyads = length(standat$combispec3_dyads)
  standat$combispec4_dyads = which(dyads$combi_num == 4)
  standat$n_combispec4_dyads = length(standat$combispec4_dyads)
  standat$combispec5_dyads = which(dyads$combi_num == 5)
  standat$n_combispec5_dyads = length(standat$combispec5_dyads)
  standat$combispec6_dyads = which(dyads$combi_num == 6)
  standat$n_combispec6_dyads = length(standat$combispec6_dyads)

  list(standat = standat, id_table = ids, dyad_table = dyads)
}


make_assoc_standat <- function(xdata, groupyear) {
  temp <- xdata[xdata$groupyear == groupyear, ]
  ids <- unique(c(temp$ind1, temp$ind2))
  n <- length(ids)
  m <- matrix(ncol = n, nrow = n, 0)
  colnames(m) <- rownames(m) <- ids
  o <- m
  i=1
  for (i in 1:nrow(temp)) {
    m[temp$ind1[i], temp$ind2[i]] <- temp$association2[i]
    o[temp$ind1[i], temp$ind2[i]] <- temp$association.total[i]
  }

  m <- m + t(m)
  m[lower.tri(m)] <- 0
  o <- o + t(o)
  m[o == 0] <- NA
  o[o == 0] <- NA
  o[lower.tri(o)] <- 0
  diag(o) <- 0
  diag(m) <- 0
  m[lower.tri(m)] <- 0

  make_stan_data_from_matrices(mats = list(asso = m),
                               obseff = list(o),
                               behav_types = "prop")
}

prep_dat_hum <- function(xdata) {

  if (interactive()) {
    message("convert HUMAN DATA table format in to list for Stan")
    message("also creates and individual-level and and a dyad-level table for post-processing")
  }

  ids <- data.frame(id = unique(c(xdata$individual1, xdata$individual2)))
  ids$num <- seq_len(nrow(ids))
  rownames(ids) <- ids$id
  ids$sex <- NA
  for (i in ids$id) {
    aux <- xdata[xdata$individual1 == i, c("individual1", "sex1")]
    if (nrow(aux) == 0) aux <- xdata[xdata$individual2 == i, c("individual2", "sex2")]
    if (aux[1, 2] == "m") ids[i, "sex"] <- "m"
    if (aux[1, 2] == "f") ids[i, "sex"] <- "f"
  }

  xdata$i1 <- t(ids[, "num", drop = FALSE])[1, xdata$individual1]
  xdata$i2 <- t(ids[, "num", drop = FALSE])[1, xdata$individual2]

  xdata$dyad <- apply(xdata[, c("i1", "i2")], 1, function(x)paste(sort(x), collapse = "_"))

  dyads <- data.frame(dyad = unique(xdata$dyad))
  dyads$num <- seq_len(nrow(dyads))
  rownames(dyads) <- dyads$dyad
  head(dyads)
  dyads$combi <- NA
  dyads$relatedness_deg <- NA
  dyads$is_spouse <- NA
  d=dyads$dyad[2]
  for (d in dyads$dyad) {
    auxids <- as.numeric(unlist(strsplit(d, "_")))
    auxsex <- sort(as.character(ids$sex[ids$num %in% auxids]))
    if (identical(auxsex, c("f", "f"))) dyads$combi[dyads$dyad == d] <- "ff"
    if (identical(auxsex, c("m", "m"))) dyads$combi[dyads$dyad == d] <- "mm"
    if (identical(auxsex, c("f", "m"))) dyads$combi[dyads$dyad == d] <- "mix"

    # kin/relatedness
    auxkin <- which((xdata$i1 == auxids[1] & xdata$i2 == auxids[2]) | (xdata$i1 == auxids[2] & xdata$i2 == auxids[1]))
    auxspouse <- xdata$partners[auxkin]
    auxkin <- xdata$kin[auxkin]
    if (length(unique(auxkin)) != 1) stop("kin error")
    if (length(unique(auxspouse)) != 1) stop("spouse error")
    dyads[d, "relatedness_deg"] <- auxkin[1]
    dyads[d, "is_spouse"] <- auxspouse[1]
  }
  table(dyads$combi, useNA = "a")
  table(dyads$relatedness_deg, useNA = "a")
  table(dyads$is_spouse, useNA = "a")

  xdata$dyad_num <- NA
  for (i in seq_len(nrow(xdata))) xdata$dyad_num[i] <- dyads$num[dyads$dyad == xdata$dyad[i]]
  all(rowSums(table(xdata$dyad_num, xdata$dyad) != 0) == 1)

  dyads$combi_num <- 1
  dyads$combi_num[dyads$combi == "mix"] <- 2
  dyads$combi_num[dyads$combi == "mm"] <- 3

  standat <- list(
    n_obs = nrow(xdata),
    n_ids = nrow(ids),
    n_dyads = nrow(dyads),
    n_behs = 4,
    id1_for_obs = xdata$i1,
    id2_for_obs = xdata$i2,
    dyad_for_obs = xdata$dyad_num,
    id_is_female = as.numeric(ids$sex == "f"),
    combi_for_dyad = dyads$combi_num,
    relatedness_deg = dyads$relatedness_deg,
    dyad_is_spouse = dyads$is_spouse,
    obseff = rep(4, nrow(xdata)),
    food = xdata$food,
    assoc = xdata$foraging_group_membership,
    help = xdata$foraging_help,
    friendship = xdata$friendship
  )

  standat$n_cors <- standat$n_behs * (standat$n_behs - 1) / 2
  standat$sex_for_id <- standat$id_is_female + 1


  standat$n_females <- sum(ids$sex == "f")
  standat$index_female_ids <- which(ids$sex == "f")
  standat$n_males <- sum(ids$sex == "m")
  standat$index_male_ids <- which(ids$sex == "m")

  standat$n_ff <- sum(dyads$combi == "ff")
  standat$index_ff_dyads <- which(dyads$combi == "ff")
  standat$n_mix <- sum(dyads$combi == "mix")
  standat$index_mix_dyads <- which(dyads$combi == "mix")
  standat$n_mm <- sum(dyads$combi == "mm")
  standat$index_mm_dyads <- which(dyads$combi == "mm")

  standat
}
