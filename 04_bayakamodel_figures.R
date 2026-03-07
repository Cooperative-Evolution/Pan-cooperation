library(cmdstanr)
pathtoresults <- "replication_pack/results/bayaka"

# reread results from disk ----
r <- readRDS(file.path(pathtoresults, "modelenv_bayaka_cooperation.rds"))
sdat <- readRDS(file.path(pathtoresults, "standat_bayaka_cooperation.rds"))
rs <- readRDS(file.path(pathtoresults, "summarydataframe_bayaka_cooperation.rds"))


# ______________________________________________________________________________
# figure 4 ---- correlations of dyad-specific propensities
# ______________________________________________________________________________

# behavior labels:
beh_labs_bayaka <- tolower(c("Sharing", "Joint Foraging", "Help", "Socializing"))
beh_combi_labs_bayaka <- tolower(c(
  paste(beh_labs_bayaka[1], "-", beh_labs_bayaka[2]),
  paste(beh_labs_bayaka[1], "-", beh_labs_bayaka[3]),
  paste(beh_labs_bayaka[2], "-", beh_labs_bayaka[3]),
  paste(beh_labs_bayaka[1], "-", beh_labs_bayaka[4]),
  paste(beh_labs_bayaka[2], "-", beh_labs_bayaka[4]),
  paste(beh_labs_bayaka[3], "-", beh_labs_bayaka[4])
))
femcol <- "#F5191C"
malecol <- "#3B99B1"
mixcol <- "gold"


x <- r$draws(format = "draws_matrix")
x <- x[, grepl("^cor", colnames(x))]
x <- t(apply(x, 2, quantile, probs = c(0.055, 0.25, 0.5, 0.75, 0.945)))
y <- data.frame(stanvar = rownames(x),
                "bcomb" = NA,
                "sex" = NA,
                check.names = FALSE)
x <- cbind(y, x)
x$bcomb <- beh_combi_labs_bayaka


n <- table(sdat$combi_for_dyad)
names(n) <- c("ff", "mix", "mm")

# go dyad
x <- x[grepl("dyad", x$stanvar), ]
x <- x[order(x$stanvar, decreasing = TRUE), ]

layout(matrix(1:4, nrow = 1), width = c(0.7, 1, 1, 1))

par(mar = c(3, 0, 1.5, 0), mgp = c(1.2, 0.3, 0), family = "serif", tcl = -0.2)
plot(0, 0, "n", xlim = c(-1, 1), ylim = c(0.5, 6.5), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
text(rep(0, 6), 1:6, labels = x$bcomb)

par(mar = c(3, 1, 1.5, 1), cex.axis = 0.9)

plot(0, 0, "n", xlim = c(-1, 1), ylim = c(0.5, 6.5), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
pdata <- x[grepl("_femfem", x$stanvar), ]
segments(0, 0.5, 0, 6.5, lty = 2, lwd = 0.5, col = "grey")
segments(x0 = pdata$`5.5%`, y0 = 1:6, x1 = pdata$`94.5%`, y1 = 1:6, lwd = 1, col = femcol)
segments(x0 = pdata$`25%`, y0 = 1:6, x1 = pdata$`75%`, y1 = 1:6, lwd = 2.5, col = femcol)
points(pdata$`50%`, 1:6, pch = 21, cex = 1.5, bg = femcol)
axis(1)
title(main = "\u2640\u2640", cex.main = 3) # female female
title(xlab = "estimated correlation coefficient")
text(1, par()$usr[4], bquote(italic(n)[dyads] == .(n["ff"])), xpd = TRUE, adj = 1, cex = 0.8)

plot(0, 0, "n", xlim = c(-1, 1), ylim = c(0.5, 6.5), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
pdata <- x[grepl("_mix", x$stanvar), ]
segments(0, 0.5, 0, 6.5, lty = 2, lwd = 0.5, col = "grey")
segments(x0 = pdata$`5.5%`, y0 = 1:6, x1 = pdata$`94.5%`, y1 = 1:6, lwd = 1, col = mixcol)
segments(x0 = pdata$`25%`, y0 = 1:6, x1 = pdata$`75%`, y1 = 1:6, lwd = 2.5, col = mixcol)
points(pdata$`50%`, 1:6, pch = 21, cex = 1.5, bg = mixcol)
axis(1)
title(main = "\u2640\u2642", cex.main = 3) # female male
title(xlab = "estimated correlation coefficient")
text(1, par()$usr[4], bquote(italic(n)[dyads] == .(n["mix"])), xpd = TRUE, adj = 1, cex = 0.8)

plot(0, 0, "n", xlim = c(-1, 1), ylim = c(0.5, 6.5), axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
pdata <- x[grepl("_malemale", x$stanvar), ]
segments(0, 0.5, 0, 6.5, lty = 2, lwd = 0.5, col = "grey")
segments(x0 = pdata$`5.5%`, y0 = 1:6, x1 = pdata$`94.5%`, y1 = 1:6, lwd = 1, col = malecol)
segments(x0 = pdata$`25%`, y0 = 1:6, x1 = pdata$`75%`, y1 = 1:6, lwd = 2.5, col = malecol)
points(pdata$`50%`, 1:6, pch = 21, cex = 1.5, bg = malecol)
axis(1)
title(main = "\u2642\u2642", cex.main = 3) # male male
title(xlab = "estimated correlation coefficient")
text(1, par()$usr[4], bquote(italic(n)[dyads] == .(n["mm"])), xpd = TRUE, adj = 1, cex = 0.8)

