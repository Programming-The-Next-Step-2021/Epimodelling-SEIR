library(Epimodelling)


context("Basic data format checking")

test_that("The correct dimensionality after model-selection with epimodellR()",{

  mod_SIR <- epimodellR("SIR")
  mod_SIS <- epimodellR("SIS")
  mod_SIRD <- epimodellR("SIRD")
  mod_SEIR <- epimodellR("SEIR")
  mod_SEIS <- epimodellR("SEIS")

  expect_true(all.equal(length(mod_SIR),
                    length(mod_SIS),
                    length(mod_SIRD),
                    length(mod_SEIR),
                    length(mod_SEIS),
                    6), TRUE)

  expect_true(all.equal(class(mod_SIR),
                        class(mod_SIS),
                        class(mod_SIRD),
                        class(mod_SEIR),
                        class(mod_SEIS),
                        "Epimodel"), TRUE)

})


test_that("The correct dimensionality after the solve()", {

  #building the model manually
  mod_SIR <- epimodellR_manual("SIR")
  mod_SIR$var_values <- c(0.8, 0.2, 0, 1)
  mod_SIR$par_values <- runif(n = 2, min = 0, max = 1)

  #solving the model
  mod_SIR_solved <- solve(mod_SIR)

  #test for class and dimensionality of the results
  expect_equal(class(mod_SIR_solved), c("Epimodelled", "matrix"))
  expect_equal(ncol(mod_SIR_solved), length(mod_SIR$variables) + 1)

})

test_that("Correct calculations based on the automatic solve()", {

  #building the model manually
  mod_SIR <- epimodellR_manual("SIR")
  mod_SIR$var_values <- c(0.8, 0.2, 0, 1)
  mod_SIR$par_values <- runif(n = 2, min = 0, max = 1)

  #solving the model
  mod_SIR_solved <- solve(mod_SIR)


  #setting up the results matrix, variables starter values and parameters
  times <- 10/0.01
  beta <- mod_SIR$par_values[1]
  gamma <- mod_SIR$par_values[2]

  res <- matrix(
    NA, nrow = times, ncol = 5, dimnames = list(NULL, c('S', 'I', 'R', "N", 'time'))
  )
  res[1, ] <- c(mod_SIR$var_values, delta_t)

  dS <- function(S, I) -beta * I * S
  dI <- function(S, I)  beta * I * S - gamma * I


  #solver
  for (i in seq(2, times)) {
    S <- res[i-1, 1]
    I <- res[i-1, 2]

    res[i, 1] <- res[i-1, 1] + delta_t * dS(S, I)
    res[i, 2] <- res[i-1, 2] + delta_t * dI(S, I)
    res[i, 5] <- delta_t * i

  }
  res[, 4] <- rep(1, times)
  res[, 3] <- 1 - res[, 1] - res[, 2]
  res

  class(res) <- c("Epimodelled", "matrix")
  expect_equal(round(mod_SIR_solved, digits = 4), round(res, digits = 4))

})







