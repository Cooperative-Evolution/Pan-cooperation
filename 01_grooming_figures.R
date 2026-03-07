library(cmdstanr)
pathtoresults <- "replication_pack/results/pan_grooming"

# reread results from disk ----
r <- readRDS(file.path(pathtoresults, "modelenv_pan_grooming.rds"))
sdat <- readRDS(file.path(pathtoresults, "standat_pan_grooming.rds"))
skeldata <- readRDS(file.path(pathtoresults, "prediction_skeleton_pan_grooming.rds"))
xdata <- readRDS(file.path(pathtoresults, "rawdata_pan_grooming.rds"))
rs <- readRDS(file.path(pathtoresults, "summarydataframe_pan_grooming.rds"))




# ______________________________________________________________________________
# figure 1 ---- grooming partner saturation for given observation effort
# ______________________________________________________________________________

## rate posteriors ----
cols <- list( bonobo = "#DB9D85", chimp = "#6CB4D9")
pdata <- r$draws("b", format = "draws_matrix")
p_fb <- density(pdata[, 1])
p_fc <- density(pdata[, 2])
p_mb <- density(pdata[, 3])
p_mc <- density(pdata[, 4])

ymax <- max(p_fb$y, p_fc$y, p_mb$y, p_mc$y)

par(family = "serif", las = 1, mgp = c(1.5, 0.5, 0), tcl = -0.2, mar = c(2.5, 2.5, 2, 1))
textcex <- 0.8
plot(0, 0, xlim = c(0, 15), ylim = c(0, ymax * 1.35), type = "n", yaxs = 'i',
     xlab = "rate parameter", ylab = "density", axes = FALSE, xaxs = 'i')

axis(1, at = seq(0, 15, by = 3))
polygon(p_fb, col = adjustcolor(cols$bonobo, 0.8))
text(2.6, 0.8, "bonobo\nfemale", cex = textcex, adj = 0)
polygon(p_mb, col = adjustcolor(cols$bonobo, 0.8))
text(1, 0.7, "bonobo\nmale", cex = textcex, adj = 1, xpd = TRUE)

polygon(p_fc, col = adjustcolor(cols$chimp, 0.8))
text(6, 0.5, "chimpanzee\nfemale", cex = textcex, adj = 0)
polygon(p_mc, col = adjustcolor(cols$chimp, 0.8))
text(10.5, 0.25, "chimpanzee\nmale", cex = textcex, adj = 0)

# ______________________________________________________________________________
## predictions for all group-years ----

pdata <- r$draws("pred_prob", format = "draws_matrix")

cols <- list(Bompusa = "#C87A8A", Kokoalongo = "#E093C3", Ekalakala = "#DB9D85",
             Kanyawara = "#98AAE1", South = "#5CBD92", East = "#38BDBB", North = "#6CB4D9")

plot(0, 0, type = "n", xlim = c(10, 35), ylim = c(0, 1), las = 1, yaxs = "i",
     xlab = "adult group size", ylab = "est. prop. group members groomed", bty = "l")

i=1
for (i in 1:nrow(skeldata)) {
  p <- c(pdata[, i])
  # col <- adjustcolor(nd$col[i], 0.5)
  col <- adjustcolor(cols[skeldata$group[i]], 0.5)
  p <- density(p)
  newy <- p$x
  newx <- p$y/5 + skeldata$groupsize[i]
  newx2 <- skeldata$groupsize[i] - p$y/5
  if (skeldata$is_female[i]) {
    polygon(newx, newy, col = col, border = NA)
    polygon(newx, newy, col = "black", border = NA, density = 50, angle = -60, lwd = 0.5)
  } else {
    polygon(newx, newy, col = col, border = NA)
    polygon(newx, newy, col = "black", border = NA, density = 50, angle = 60, lwd = 0.5)
  }
}

for (i in 1:nrow(skeldata)) {
  med <- median(c(pdata[, i]))
  pch <- c(21, 23)[(skeldata$species[i] == "chimp") + 1]
  col <- ifelse(skeldata$is_female[i] == 1, "darkred", "darkblue")
  points(skeldata$groupsize[i], med, pch = pch, col = col, bg = cols[[skeldata$group[i]]], lwd = 4, cex = 1.5)
}
# box()

text(24.5, 0.4, "male", col = "darkblue", font = 2, adj = 0)
segments(x0 = 24.3, x1 = 22.5, y0 = 0.4, y1 = 0.43)

text(27.5, 0.1, "female", col = "darkred", font = 2, adj = 1)
segments(x0 = 27.7, x1 = 28.5, y0 = 0.1, y1 = 0.08)

### legend for groups... ----
bgroups <- c("Ekalakala", "Kokoalongo", "Bompusa")
cgroups <- c("North", "East", "South", "Kanyawara")

x <- matrix(c(bgroups, NA, cgroups), ncol = 4, byrow = TRUE)
x <- rbind(NA, x[1, ], NA, x[2, ])
x[1, 1] <- "bonobo"
x[3, 1] <- "chimpanzee"


colmat <- matrix(unlist(c(cols[bgroups], NA, cols[cgroups])), ncol = 4, byrow = TRUE)
colmat <- rbind(NA, colmat[1, ], NA, colmat[2, ])
x[x == "Bompusa"] <- "Bompusa West"

pchmat <- matrix(c(NA, NA, NA, NA,
                   21, 21, 21, NA,
                   NA, NA, NA, NA,
                   23, 23, 23, 23), ncol = 4, byrow = TRUE
)
fmat <- matrix(ncol = 4, nrow = 4, 1)
fmat[1, 1] <- fmat[3, 1] <- 2
legend("topright",
  legend = x, text.font = fmat, xpd = TRUE,
  bty = "n", ncol = 4, cex = 0.7, text.width = 2.5,
  pch = pchmat, pt.bg = colmat, pt.cex = 1.3, pt.lwd = 2, y.intersp = 1.2)

