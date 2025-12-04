test_that("mu is returned correctly", {
  result <- myncurve(1, 2, 3)
  expect_equal(result$mu, 1)
})
test_that("sigma is returned correctly", {
  result <- myncurve(1, 2, 3)
  expect_equal(result$sigma, 2)
})
test_that("probability is correctly calculated", {
  result <- myncurve(1, 2, 3)
  expected <- round(pnorm(3, mean=1, sd=2), 4)
  expect_equal(result$prob, expected)
})
