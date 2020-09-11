context("test-haystack_sparse")

# convert data to sparse matrices
dat.expression.sparse <- as(dat.expression, "CsparseMatrix")
dat.detection.sparse <- dat.expression.sparse > 1

# Test haystack 2D
set.seed(123)
x <- haystack(dat.tsne, detection = dat.detection.sparse, method = "2D")

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), "results")
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 4))
  expect_equal(x$results["gene_1", "D_KL"], 0.06756395, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "log.p.vals"], -0.003983411, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "T.counts"], 27)
})

x <- haystack(dat.tsne, detection = dat.detection.sparse, method = "highD", grid.method = "centroid", grid.points = 60)

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "grid.coordinates"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 4))
  expect_equal(class(x$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$grid.coordinates), c(60, 2))
  expect_equal(x$results["gene_1", "D_KL"], 0.6742179, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "log.p.vals"], -0.06102772, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "T.counts"], 27)
})

x <- haystack(dat.tsne, detection = dat.detection.sparse, method = "highD", grid.method = "seeding", grid.points = 60)

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "grid.coordinates"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 4))
  expect_equal(class(x$grid.coordinates)[1], "matrix")
  expect_equal(dim(x$grid.coordinates), c(60, 2))
  expect_equal(x$results["gene_1", "D_KL"], 0.4650268, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "log.p.vals"], -0.01483579, tolerance = 1e-6)
  expect_equal(x$results["gene_1", "T.counts"], 27)
})

