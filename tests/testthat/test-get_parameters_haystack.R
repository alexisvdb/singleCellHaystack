context("test-get_parameters_haystack")

y <- dat.tsne
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])

test_that("get_parameters_haystack returns valid object", {
  expect_type(param, "list")
  expect_equal(names(param), c("grid.points", "limits", "bandwidths", "dens.x", "dens.y"))
  expect_equal(param$grid.points, c(43, 41))
})
