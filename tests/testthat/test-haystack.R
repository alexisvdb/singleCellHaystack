context("test-haystack")

# Test haystack 2D
set.seed(123)
x <- haystack(dat.tsne, detection = dat.expression > 1, method = "2D")

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), "results")
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(x$results["gene_1", "D_KL"], 0.06756395)
  expect_equal(x$results["gene_1", "log.p.vals"], 0)
  expect_equal(x$results["gene_1", "T.counts"], 27)
})

x <- haystack(dat.tsne, detection = dat.expression > 1, method = "highD", grid.method = "centroid")

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "grid.coordinates"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$grid.coordinates), "matrix")
  expect_equal(dim(x$grid.coordinates), c(50, 2))
  expect_equal(x$results["gene_1", "D_KL"], 0.01781055)
  expect_equal(x$results["gene_1", "log.p.vals"], 0)
  expect_equal(x$results["gene_1", "T.counts"], 27)

})

x <- haystack(dat.tsne, detection = dat.expression > 1, method = "highD", grid.method = "seeding")

test_that("haystack works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "grid.coordinates"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(500, 3))
  expect_equal(class(x$grid.coordinates), "matrix")
  expect_equal(dim(x$grid.coordinates), c(50, 2))
  expect_equal(x$results["gene_1", "D_KL"], 0.01789938)
  expect_equal(x$results["gene_1", "log.p.vals"], 0)
  expect_equal(x$results["gene_1", "T.counts"], 27)
})
