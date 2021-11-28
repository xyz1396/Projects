context("peak_calculator_DIY")

test_that("peak_calculator_DIY", {
  print(peak_calculator_DIY("SRKSD", "N15", 0.5))
  expect_length(peak_calculator_DIY("SRKSD", "N15", 0.5), 2)
})