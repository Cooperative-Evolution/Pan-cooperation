library(cmdstanr)
pathtoresults <- "replication_pack/results/pan_associations/"


# reread results from disk ----
resmat_intercept <- readRDS(file.path(pathtoresults, "postdraws_baseline.rds"))
resmat_indi <- readRDS(file.path(pathtoresults, "postdraws_indivar.rds"))
resmat_dyad <- readRDS(file.path(pathtoresults, "postdraws_dyadvar.rds"))

# ______________________________________________________________________________
# figure 3 ---- association propensities
# ______________________________________________________________________________


cols <- c(Bompusa = "#C87A8A", Kokoalongo = "#E093C3", Ekalakala = "#DB9D85",
          Kanyawara = "#98AAE1", South = "#5CBD92", East = "#38BDBB", North = "#6CB4D9")


bonobos <- c("Bompusa", "Kokoalongo", "Ekalakala")

# layout
lmat <- matrix(c(4, 4, 4, 1, 2, 3), ncol = 3, byrow = TRUE)
layout(lmat, heights = c(1, 6))


par(mar = c(2, 4, 1, 1), mgp = c(1.7, 0.4, 0), tcl = -0.2, family = "serif")

hspread <- 28
plot(0, 0, ylim = c(0, 1), xlim = c(-hspread, hspread), yaxs = "i", type = "n",
     ylab = "estimated dyadic baseline association", xlab = "", las = 1, bty = "l",
     xaxs = "i", axes = FALSE)
title(xlab = "density", line = 0.5)
axis(2, las = 1)
box(bty = "l")
abline(v = 0)
text(-hspread, par()$usr[4], "(A)", adj = c(-0.5, 0.5), xpd = TRUE)

i=colnames(resmat_intercept)[1]
for (i in colnames(resmat_intercept)) {
  g <- strsplit(i, " ")[[1]][1]
  # flip
  pp <- plogis(resmat_intercept[, i])
  p <- density(pp)
  x <- p$y
  y <- p$x
  p$y <- y
  p$x <- x
  if (g %in% bonobos) p$x <- p$x * (-1)

  col <- adjustcolor(cols[g], 0.7)
  polygon(p, col = col, border = grey(0.2), lwd = 0.5)
}
abline(v = 0)
text(x = c(-0.5 * hspread, 0.5 * hspread), y = par()$usr[4] * 0.92, labels = c("bonobo", "chimpanzee"), cex = 1.1, adj = c(0.5, 0))

for (i in colnames(resmat_intercept)) {
  pp <- plogis(resmat_intercept[, i])
  multipl <- 1
  if (strsplit(i, " ")[[1]][1] %in% bonobos) multipl <- -1
  segments(x0 = 0, y0 = median(pp), x1 = hspread/10 * multipl, y1 = median(pp), lwd = 1.5)
}

## dyad propensities ----
hspread <- 30
plot(0, 0, ylim = c(0, 1.6), xlim = c(-hspread, hspread), yaxs = "i", type = "n",
     ylab = "estimated variation in\ndyadic association propensities", xlab = "",
     las = 1, bty = "l", xaxs = "i", axes = FALSE)
title(xlab = "density", line = 0.5)
axis(2, las = 1)
box(bty = "l")
abline(v = 0)
text(-hspread, par()$usr[4], "(B)", adj = c(-0.5, 0.5), xpd = TRUE)

i=colnames(resmat_dyad)[1]
for (i in colnames(resmat_dyad)) {
  g <- strsplit(i, " ")[[1]][1]
  # flip
  pp <- resmat_dyad[, i]
  p <- density(pp)
  x <- p$y
  y <- p$x
  p$y <- y
  p$x <- x
  if (g %in% bonobos) p$x <- p$x * (-1)
  col <- adjustcolor(cols[g], 0.7)
  polygon(p, col = col, border = grey(0.2), lwd = 0.5)
}
abline(v = 0)
text(x = c(-0.5 * hspread, 0.5 * hspread), y = par()$usr[4] * 0.92, labels = c("bonobo", "chimpanzee"), cex = 1.1, adj = c(0.5, 0))

for (i in colnames(resmat_dyad)) {
  pp <- resmat_dyad[, i]
  multipl <- 1
  if (strsplit(i, " ")[[1]][1] %in% bonobos) multipl <- -1
  segments(x0 = 0, y0 = median(pp), x1 = hspread/10 * multipl, y1 = median(pp), lwd = 1.5)

}



## indi propensities ----
hspread <- 13
plot(0, 0, ylim = c(0, 1.6), xlim = c(-hspread, hspread), yaxs = "i", type = "n",
     ylab = "estimated variation in\nindividual association propensities", xlab = "",
     las = 1, bty = "l", xaxs = "i", axes = FALSE, xpd = TRUE)
title(xlab = "density", line = 0.5)
axis(2, las = 1)
box(bty = "l")
abline(v = 0)
text(-hspread, par()$usr[4], "(C)", adj = c(-0.5, 0.5), xpd = TRUE)

i=colnames(resmat_indi)[1]
for (i in colnames(resmat_indi)) {
  g <- strsplit(i, " ")[[1]][1]
  # flip
  pp <- resmat_indi[, i]
  p <- density(pp)
  x <- p$y
  y <- p$x
  p$y <- y
  p$x <- x
  if (g %in% bonobos) {
    p$x <- p$x * (-1)
  }
  col <- adjustcolor(cols[g], 0.7)
  polygon(p, col = col, border = grey(0.2), lwd = 0.5)
}
abline(v = 0)
text(x = c(-0.5 * hspread, 0.5 * hspread), y = par()$usr[4] * 0.92, labels = c("bonobo", "chimpanzee"), cex = 1.1, adj = c(0.5, 0))

for (i in colnames(resmat_indi)) {
  pp <- resmat_indi[, i]
  multipl <- 1
  if (strsplit(i, " ")[[1]][1] %in% bonobos) multipl <- -1
  segments(x0 = 0, y0 = median(pp), x1 = hspread/10 * multipl, y1 = median(pp), lwd = 1.5)

}

# legend(s)----
par(mar = c(0.5, 0, 0, 0))
plot(0, 0, axes = F, ann = F, type = "n")

legend(x = -0.3, y = 1, xjust = 0.5,
       legend = c("Bompusa West", "Kokoalongo", "Ekalakala", NA),
       ncol = 2, cex = 0.9, pch = 15, col = adjustcolor(c(cols[c("Bompusa", "Kokoalongo", "Ekalakala")], NA), 1),
       pt.cex = 1.3, xpd = TRUE, title = expression(italic("P. paniscus")),
       bty = 'o', , bg = grey(0.9), box.col = NA
)

legend(x = 0.3, y = 1, xjust = 0.5,
       legend = c("Kanyawara", "South", "East", "North"),
       ncol = 2, cex = 0.9, pch = 15, col = adjustcolor(c(cols[c("Kanyawara", "South", "East", "North")]), 1),
       pt.cex = 1.3, xpd = TRUE, title = expression(italic("P. troglodytes")),
       bty = 'o', bg = grey(0.9), box.col = NA
)

