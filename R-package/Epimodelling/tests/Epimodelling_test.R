library(Epimodelling)


context("Epimodelling basic functionality")

test_that("check the class and length of the list in Epimodeller()", {

  t_epimodeller <- epimodellR("SIR")

  expect_equal(length(t_epimodeller), 6)
  expect_equal(class(t_epimodeller), "Epimodel")

})

#expand this test more
