context("peak_calculator")

test_that("peak_calculator works", {
  print(peak_calculator("SRKSD"))
  expect_length(peak_calculator("SRKSD"), 2)
})