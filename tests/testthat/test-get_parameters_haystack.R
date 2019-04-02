context("test-get_parameters_haystack")

set.seed(123)
y <- Rtsne::Rtsne(mtcars, perplexity = 10)$Y
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])

test_that("get_parameters_haystack returns valid object", {
  expect_type(param, "list")
  expect_equal(names(param), c("grid.points", "limits", "bandwidths", "dens.x", "dens.y"))
  expect_equal(param$grid.points, c(40, 51))
})
