
plot_post_cor_4axes <- function(model_env,
                                standat,
                                aggr = FALSE,
                                by_individual = TRUE,
                                add_pplus = FALSE
) {
  if (interactive()) {
    message("posteriors of the correlation coefficients")
  }

  pairlabs <- matrix(c("grooming", "aggression", # vec[1]
                       "grooming", "coalition", # vec[2]
                       "aggression", "coalition", # vec [3]
                       "grooming", "sharing", # vec[4]
                       "aggression", "sharing", # vec[5]
                       "coalition", "sharing" # vec[6]
  ), ncol = 2, byrow = TRUE)
  pairlabs <- apply(pairlabs, 1, paste, collapse = "-")


  if (by_individual) {
    # 1 = bonobos
    d_f1 <- data.frame(model_env$summary("cors_indi_sexspec1", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_m1 <- data.frame(model_env$summary("cors_indi_sexspec2", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    # 2 = chimps
    d_f2 <- data.frame(model_env$summary("cors_indi_sexspec3", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_m2 <- data.frame(model_env$summary("cors_indi_sexspec4", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))

    # pplus
    pplus_fem <- sapply(1:6, function(i) {
      aux <- model_env$draws(c("cors_indi_sexspec1", "cors_indi_sexspec3"), format = "draws_matrix")
      aux <- aux[, grepl(paste0("[", i, "]"), colnames(aux), fixed = TRUE)]
      mean(c(aux[, 2] - aux[, 1]) > 0)
    })
    pplus_male <- sapply(1:6, function(i) {
      aux <- model_env$draws(c("cors_indi_sexspec2", "cors_indi_sexspec4"), format = "draws_matrix")
      aux <- aux[, grepl(paste0("[", i, "]"), colnames(aux), fixed = TRUE)]
      mean(c(aux[, 2] - aux[, 1]) > 0)
    })

    names(pplus_fem) <- pairlabs
    names(pplus_male) <- pairlabs
    d_f1$variable <- pairlabs
    d_m1$variable <- pairlabs
    d_f2$variable <- pairlabs
    d_m2$variable <- pairlabs

    # split by aggression
    if (!is.null(aggr)) {
      if (aggr) {
        d_f1 <- d_f1[grepl("aggression", d_f1$variable), ]
        d_m1 <- d_m1[grepl("aggression", d_m1$variable), ]
        d_f2 <- d_f2[grepl("aggression", d_f2$variable), ]
        d_m2 <- d_m2[grepl("aggression", d_m2$variable), ]
        pplus_fem <- pplus_fem[grepl("aggression", names(pplus_fem))]
        pplus_male <- pplus_male[grepl("aggression", names(pplus_male))]
      } else {
        d_f1 <- d_f1[!grepl("aggression", d_f1$variable), ]
        d_m1 <- d_m1[!grepl("aggression", d_m1$variable), ]
        d_f2 <- d_f2[!grepl("aggression", d_f2$variable), ]
        d_m2 <- d_m2[!grepl("aggression", d_m2$variable), ]
        pplus_fem <- pplus_fem[!grepl("aggression", names(pplus_fem))]
        pplus_male <- pplus_male[!grepl("aggression", names(pplus_male))]

      }
    }
    outres <- list(bon_fem = d_f1, bon_male = d_m1, chi_fem = d_f2, chi_male = d_m2)

    n_cors_tot <- nrow(d_f1) * 4
    n_cors_per <- nrow(d_f1)

    # color defs ("zissou" palette)
    xcols <- c(rep("#F5191C", n_cors_per * 2), rep("#3B99B1", n_cors_per * 2))

    op <- par(no.readonly = TRUE)
    lmat <- matrix(c(8, 6, 7, 5, 1, 3, 5, 2, 4), ncol = 3)
    layout(lmat, widths = c(1, 4, 4), heights = c(1, 4, 4))
    i=2
    for (i in 1:4) {
      if (i == 2 || i == 4) {
        par(mar = c(3, 3, 1, 4), mgp = c(1.2, 0.35, 0), tcl = -0.2)
      } else {
        par(mar = c(3, 7, 1, 0), mgp = c(1.2, 0.35, 0), tcl = -0.2)
      }
      if (i == 1) {
        dd <- d_f1
        ni <- standat$n_sexspec1_ids
      }
      if (i == 2) {
        dd <- d_f2
        ni <- standat$n_sexspec3_ids
        dd$pplus <- sprintf("%0.2f", pplus_fem)
      }
      if (i == 3) {
        dd <- d_m1
        ni <- standat$n_sexspec2_ids
      }
      if (i == 4) {
        dd <- d_m2
        ni <- standat$n_sexspec4_ids
        dd$pplus <- sprintf("%0.2f", pplus_male)
      }
      # reverse order in data frame so that groom_coa is at the top
      dd <- dd[rev(1:nrow(dd)), ]
      xcol <- xcols[seq(0, n_cors_tot + 1, by = n_cors_per)[i] + (1:n_cors_per)]
      plot(0, 0, xlim = c(-1, 1), ylim = c(0, n_cors_per + 1), xaxs = "i", yaxs = "i",
           type = "n", axes = FALSE, xlab = "estimated correlation coefficient", ylab = "", cex.lab = 0.9)
      abline(v = 0, lty = 2, col = "grey", lwd = 0.5)
      segments(dd$X5.5., 1:n_cors_per, dd$X94.5., 1:n_cors_per, col = xcol)
      segments(dd$X25., 1:n_cors_per, dd$X75., 1:n_cors_per, lwd = 3, col = xcol)
      points(dd$X50., 1:n_cors_per, pch = 21, cex = 2, bg = xcol)
      if (!(i == 2 || i == 4)) {
        axis(2, at = 1:n_cors_per, labels = dd$variable,
             las = 1, tcl = 0, lwd = 0)
      } else {
        # add pplus
        for (kk in 1:3) {
          if (add_pplus) {
            text(1, kk, substitute(italic(p)^"+" == value,
                                   list(value = dd$pplus[kk])),
                 xpd = TRUE, adj = 0, cex = 0.8)
          }
        }
      }

      axis(1)
      usr <- par("usr")
      legend_x <- usr[2] - 0.35 * (usr[2] - usr[1])
      legend_y <- usr[4] + 0.1 * (usr[4] - usr[3])
      legend(x = legend_x, y = legend_y, pch = NA, legend = paste("n =", ni, "indiv."),
             cex = 0.8, bty = "n", xpd = TRUE, horiz = TRUE)

    }
    par(mar = c(0, 3.2, 0, 1))
    plot(0, 0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
    text(-0.5, -0.75, "bonobo", cex = 1.5, adj = c(0.5, 0), xpd = TRUE)
    text(0.5, -0.75, "chimpanzee", cex = 1.5, adj = c(0.5, 0))
    yloc <- 0
    xxcol <- "#4D4D4D" # grey(0.3)
    points(0, yloc, pch = 16, cex = 2, col = xxcol)
    fac <- 0.5
    segments(x0 = -0.25 * fac, y0 = yloc, x1 = 0.25 * fac, y1 = yloc, lwd = 2.6, col = xxcol)
    segments(x0 = -0.445 * fac, y0 = yloc, x1 = 0.445 * fac, y1 = yloc, lwd = 1, col = xxcol)
    text(0, -0.5, "median", adj = c(0.5, 1), cex = 0.8, col = xxcol)
    text(fac/5, -0.5, "50%", adj = c(0.5, 1), cex = 0.8, col = xxcol)
    text(fac/2.5, -0.5, "89%", adj = c(0.5, 1), cex = 0.8, col = xxcol)

    plot(0, 0, "n", axes = FALSE, ann = FALSE)
    text(-0.2, 0.2, "\u2640", cex = 3.5, font = 2)
    plot(0, 0, "n", axes = FALSE, ann = FALSE)
    text(-0.2, 0.2, "\u2642", cex = 3.5, font = 2)

  }

  if (!by_individual) {
    # bonobos
    d_ff1 <- data.frame(model_env$summary("cors_dyad_combispec1", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_mix1 <- data.frame(model_env$summary("cors_dyad_combispec2", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_mm1 <- data.frame(model_env$summary("cors_dyad_combispec3", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    # chimps
    d_ff2 <- data.frame(model_env$summary("cors_dyad_combispec4", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_mix2 <- data.frame(model_env$summary("cors_dyad_combispec5", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))
    d_mm2 <- data.frame(model_env$summary("cors_dyad_combispec6", ~quantile(.x, probs = c(0.055, 0.25, 0.5, 0.75, 0.945))))

    d_ff1$variable <- pairlabs
    d_mix1$variable <- pairlabs
    d_mm1$variable <- pairlabs
    d_ff2$variable <- pairlabs
    d_mix2$variable <- pairlabs
    d_mm2$variable <- pairlabs

    # pplus
    pplus_ff <- sapply(1:6, function(i) {
      aux <- model_env$draws(c("cors_dyad_combispec1", "cors_dyad_combispec4"), format = "draws_matrix")
      aux <- aux[, grepl(paste0("[", i, "]"), colnames(aux), fixed = TRUE)]
      mean(c(aux[, 2] - aux[, 1]) > 0)
    })
    pplus_mix <- sapply(1:6, function(i) {
      aux <- model_env$draws(c("cors_dyad_combispec2", "cors_dyad_combispec5"), format = "draws_matrix")
      aux <- aux[, grepl(paste0("[", i, "]"), colnames(aux), fixed = TRUE)]
      mean(c(aux[, 2] - aux[, 1]) > 0)
    })
    pplus_mm <- sapply(1:6, function(i) {
      aux <- model_env$draws(c("cors_dyad_combispec3", "cors_dyad_combispec6"), format = "draws_matrix")
      aux <- aux[, grepl(paste0("[", i, "]"), colnames(aux), fixed = TRUE)]
      mean(c(aux[, 2] - aux[, 1]) > 0)
    })

    names(pplus_ff) <- pairlabs
    names(pplus_mix) <- pairlabs
    names(pplus_mm) <- pairlabs

    # split by aggression
    if (!is.null(aggr)) {
      if (aggr) {
        d_ff1 <- d_ff1[grepl("aggression", d_ff1$variable), ]
        d_mix1 <- d_mix1[grepl("aggression", d_mix1$variable), ]
        d_mm1 <- d_mm1[grepl("aggression", d_mm1$variable), ]
        d_ff2 <- d_ff2[grepl("aggression", d_ff2$variable), ]
        d_mix2 <- d_mix2[grepl("aggression", d_mix2$variable), ]
        d_mm2 <- d_mm2[grepl("aggression", d_mm2$variable), ]
        pplus_ff <- pplus_ff[grepl("aggression", names(pplus_ff))]
        pplus_mix <- pplus_mix[grepl("aggression", names(pplus_mix))]
        pplus_mm <- pplus_mm[grepl("aggression", names(pplus_mm))]

      } else {
        d_ff1 <- d_ff1[!grepl("aggression", d_ff1$variable), ]
        d_ff1 <- d_ff1[!grepl("aggression", d_ff1$variable), ]
        d_mix1 <- d_mix1[!grepl("aggression", d_mix1$variable), ]
        d_mm1 <- d_mm1[!grepl("aggression", d_mm1$variable), ]
        d_ff2 <- d_ff2[!grepl("aggression", d_ff2$variable), ]
        d_mix2 <- d_mix2[!grepl("aggression", d_mix2$variable), ]
        d_mm2 <- d_mm2[!grepl("aggression", d_mm2$variable), ]
        pplus_ff <- pplus_ff[!grepl("aggression", names(pplus_ff))]
        pplus_mix <- pplus_mix[!grepl("aggression", names(pplus_mix))]
        pplus_mm <- pplus_mm[!grepl("aggression", names(pplus_mm))]
      }
    }
    outres <- list(bon_ff = d_ff1, bon_mix = d_mix1, bon_mm = d_mm1,
                   chi_ff = d_ff2, chi_mix = d_mix2, chi_mm = d_mm2)

    n_cors_tot <- nrow(d_ff1) * 6
    n_cors_per <- nrow(d_ff1)

    xcols <- c(rep("#F5191C", n_cors_per * 2), rep("gold", n_cors_per * 2), rep("#3B99B1", n_cors_per * 2))

    op <- par(no.readonly = TRUE)
    lmat <- matrix(c(11, 8, 9, 10, 7, 1, 3, 5, 7, 2, 4, 6), ncol = 3)
    layout(lmat, widths = c(1.9, 4, 4), heights = c(1, 4, 4, 4))
    for (i in 1:6) { # six panels
      # Set margins based on the value of i
      if (i == 2 || i == 4 || i == 6) {
        par(mar = c(3, 3, 1, 4), mgp = c(1.2, 0.35, 0), tcl = -0.2)
      } else {
        par(mar = c(3, 7, 1, 0), mgp = c(1.2, 0.35, 0), tcl = -0.2)
      }
      if (i == 1) {
        dd <- d_ff1
        ni <- standat$n_combispec1_dyads
      }
      if (i == 2) {
        dd <- d_ff2
        ni <- standat$n_combispec4_dyads
        dd$pplus <- sprintf("%0.2f", pplus_ff)
      }
      if (i == 3) {
        dd <- d_mix1
        ni <- standat$n_combispec2_dyads
      }
      if (i == 4) {
        dd <- d_mix2
        ni <- standat$n_combispec5_dyads
        dd$pplus <- sprintf("%0.2f", pplus_mix)
      }
      if (i == 5) {
        dd <- d_mm1
        ni <- standat$n_combispec3_dyads
      }
      if (i == 6) {
        dd <- d_mm2
        ni <- standat$n_combispec6_dyads
        dd$pplus <- sprintf("%0.2f", pplus_mm)
      }

      # reverse order in data frame so that groom_coa is at the top
      dd <- dd[rev(1:nrow(dd)), ]

      xcol <- xcols[seq(0, n_cors_tot + 1, by = n_cors_per)[i] + (1:n_cors_per)]
      plot(0, 0, xlim = c(-1, 1), ylim = c(0, n_cors_per + 1), xaxs = "i", yaxs = "i",
           type = "n", axes = FALSE, xlab = "estimated correlation coefficient", ylab = "", cex.lab = 0.9)
      abline(v = 0, lty = 2, col = "grey", lwd = 0.5)
      segments(dd$X5.5., 1:n_cors_per, dd$X94.5., 1:n_cors_per, col = xcol)
      segments(dd$X25., 1:n_cors_per, dd$X75., 1:n_cors_per, lwd = 3, col = xcol)
      points(dd$X50., 1:n_cors_per, pch = 21, cex = 2, bg = xcol)
      axis(1)
      if (!(i == 2 || i == 4 || i == 6)) {
        axis(2, at = 1:n_cors_per, labels = dd$variable,
             las = 1, tcl = 0, lwd = 0)
      } else {
        # add pplus
        for (kk in 1:3) {
          if (add_pplus) {
            text(1, kk, substitute(italic(p)^"+" == value,
                                   list(value = dd$pplus[kk])),
                 xpd = TRUE, adj = 0, cex = 0.8)
          }

        }
      }
      usr <- par("usr")
      legend_x <- usr[2] - 0.35 * (usr[2] - usr[1])
      legend_y <- usr[4] + 0.12 * (usr[4] - usr[3])
      legend(x = legend_x, y = legend_y, pch = NA, legend = paste("n =", ni, "dyads"),
             cex = 0.8, bty = "n", xpd = TRUE, horiz = TRUE)

    }
    par(mar = c(0, 3.2, 0, 1))
    plot(0, 0, xlim = c(-1, 1), ylim = c(-1, 1), type = "n", ann = FALSE, axes = FALSE, xaxs = "i", yaxs = "i")
    text(-0.5, -0.75, "bonobo", cex = 1.5, adj = c(0.5, 0), xpd = TRUE)
    text(0.5, -0.75, "chimpanzee", cex = 1.5, adj = c(0.5, 0))
    yloc <- 0
    xxcol <- "#4D4D4D" # grey(0.3)
    points(0, yloc, pch = 16, cex = 2, col = xxcol)
    fac <- 0.5
    segments(x0 = -0.25 * fac, y0 = yloc, x1 = 0.25 * fac, y1 = yloc, lwd = 2.6, col = xxcol)
    segments(x0 = -0.445 * fac, y0 = yloc, x1 = 0.445 * fac, y1 = yloc, lwd = 1, col = xxcol)
    text(0, -0.5, "median", adj = c(0.5, 1), cex = 0.8, col = xxcol)
    text(fac/5, -0.5, "50%", adj = c(0.5, 1), cex = 0.8, col = xxcol)
    text(fac/2.5, -0.5, "89%", adj = c(0.5, 1), cex = 0.8, col = xxcol)

    plot(0, 0, "n", axes = FALSE, ann = FALSE)
    text(-0.2, 0.2, "\u2640\u2640", cex = 3.5, font = 2)
    plot(0, 0, "n", axes = FALSE, ann = FALSE)
    text(-0.2, 0.2, "\u2640\u2642", cex = 3.5, font = 2)
    plot(0, 0, "n", axes = FALSE, ann = FALSE)
    text(-0.2, 0.2, "\u2642\u2642", cex = 3.5, font = 2)
  }

  par(op)

  invisible(outres)
}

