
simple_trace <- function(mod, par) {
  xcols <- hcl.colors(4, "viridis")
  temp <- mod$draws(par)
  temp <- apply(temp, 2, \(x)x)
  plot (0, 0, type = "n", xlim = c(0, length(temp)/ncol(temp)), ylim = range(temp), axes = FALSE, ann = FALSE)
  sapply(seq_len(ncol(temp)), \(i) {
    # points(p[, 2, ][, , i], type = "l", col = hcl.colors(4, "zissou")[2], lwd = 0.5)
    points(temp[, i], type = "l", col = xcols[i], lwd = 0.5)
    invisible(NULL)
  })
  axis(2, las = 1)
  box(bty = "l")
}
