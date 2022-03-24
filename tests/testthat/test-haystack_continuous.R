context("test-continuous")

# Test haystack 2D
set.seed(123)
x <- haystack_continuous_2D(x = dat.tsne[,1], y = dat.tsne[,2], expression = dat.expression)

test_that("haystack continuous 2D works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(x$results["gene_24", "D_KL"], 1.74568845, tolerance = 1e-6)
  expect_equal(x$results["gene_24", "log.p.vals"], -7.79657, tolerance = 1e-6)
})

x <- haystack_continuous_highD(x = dat.tsne, expression = dat.expression, grid.points = 60)

test_that("haystack continuous highD works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$info$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$info$grid.coordinates), c(60, 2))
  expect_equal(x$results["gene_24", "D_KL"], 2.459008, tolerance = 1e-6)
  expect_equal(x$results["gene_24", "log.p.vals"], -28.16831, tolerance = 1e-6)
})

x <- haystack_continuous_highD(x = dat.tsne, expression = dat.expression, grid.method = "seeding", grid.points = 60)

test_that("haystack continuous highD works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$info$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$info$grid.coordinates), c(60, 2))
  expect_equal(x$results["gene_24", "D_KL"], 2.448916, tolerance = 1e-6)
  expect_equal(x$results["gene_24", "log.p.vals"], -22.1835, tolerance = 1e-6)
})
