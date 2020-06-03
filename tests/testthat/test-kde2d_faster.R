context("test-kde2d_faster")

y <- dat.tsne
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])
z <- singleCellHaystack:::kde2d_faster(param$dens.x, param$dens.y)

test_that("kde2d_faster works", {
  expect_type(z, "double")
  expect_equal(class(z)[1], "matrix")
  expect_equal(dim(z), c(43, 41))
  expect_equal(dim(z), param$grid.points)
})
