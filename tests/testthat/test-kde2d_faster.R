context("test-kde2d_faster")

set.seed(123)
y <- Rtsne::Rtsne(mtcars, perplexity = 10)$Y
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])
z <- singleCellHaystack:::kde2d_faster(param$dens.x, param$dens.y)

test_that("kde2d_faster works", {
  expect_type(z, "double")
  expect_equal(class(z)[1], "matrix")
  expect_equal(dim(z), c(40, 51))
  expect_equal(dim(z), param$grid.points)
})
