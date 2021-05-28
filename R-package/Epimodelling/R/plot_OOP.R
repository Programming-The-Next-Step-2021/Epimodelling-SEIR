#'Method for plotting Epimodels
#'
#'The method plots the results of the solved system of differential equations based on the \emph{Epimodel} object.
#'
#'@param res The results of the \emph{solve()} method, which is an \emph{Epimodelled} object
#'@param main Title for the plot
#'
#'@examples
#'
#'mod <- epimodellR("SIR")
#'mod_solved <- solve(mod)
#'plot(mod_solved)
#'
#' @export

plot.Epimodelled <- function(res, main = ""){

  vars <- colnames(res)
  vars <- vars[-(length(vars) : (length(vars) - 1)) ]

  time <- res[, ncol(res)]

  matplot(
    time,
    res[, 1:(ncol(res) - 2)],
    type = 'l',
    axes = FALSE,
    lty = 1, lwd = 2,
    ylab = 'Population %', xlab = 'Days', ylim = c(0, 1), xlim = c(0, tail(time, 1)),
    main = main, cex.main = 1.75, cex.lab = 1.25, font.main = 1,
    xaxs = 'i', yaxs = 'i'
  )

  axis(1, cex.axis = 1.5)
  axis(2, las = 2, cex.axis = 1.5)

  legend("topright",
         legend = vars,
         lty = 1, lwd = 2,
         bty = 'n', cex = 1.5,
         col = 1:length(vars)
  )


}


