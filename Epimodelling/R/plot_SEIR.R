#' @export

plot_SIRS <- function(res, main = '') {
  time <- res[, 6]

  matplot(
    time, res[, 2:5]/res[1,1], type = 'l', axes = FALSE, lty = 1, lwd = 2,
    ylab = 'Subpopulations', xlab = 'Days', ylim = c(0, 1),
    main = main, cex.main = 1.75, cex.lab = 1.25,
    font.main = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, tail(time, 1))
  )

  axis(1, cex.axis = 1.5)
  axis(2, las = 2, cex.axis = 1.5)
  legend(
    30, 0.95, legend = c('S', 'I', 'R'),
    lty = 1, lwd = 2, bty = 'n', cex = 1.5
  )
}


