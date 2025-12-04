context("test-haystack_sparse")

# convert data to sparse matrices
dat.expression.sparse <- as(dat.expression, "CsparseMatrix")

set.seed(123)

suppressWarnings({
x <- haystack(dat.tsne, dat.expression.sparse, grid.method = "centroid", grid.points = 60)
})

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$info$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$info$grid.coordinates), c(60, 2))
  #expect_equal(x$results["gene_1", "D_KL"], 0.3401927, tolerance = 1e-6)
  #expect_equal(x$results["gene_1", "log.p.vals"], -0.4574128, tolerance = 1e-6)
})

suppressWarnings({
x <- haystack(dat.tsne, dat.expression.sparse, grid.method = "seeding", grid.points = 60)
})

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$info$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$info$grid.coordinates), c(60, 2))
  #expect_equal(x$results["gene_1", "D_KL"], 0.2716761, tolerance = 1e-6)
  #expect_equal(x$results["gene_1", "log.p.vals"], -0.4350106, tolerance = 1e-6)
})

