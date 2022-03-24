context("test-spatial")

# Test haystack spatial
set.seed(123)
# quickly make a small artificial dataset
spots <- 500
genes <- 100
coordinates <- matrix(runif(n=spots*2), ncol=2)
counts <- matrix(rpois(n=spots*genes, lambda = 1.5), nrow=genes)
rownames(counts) <- paste0("gene_",1:genes)
# add 1 gene that has high expression around the 0.5 / 0.5 point
center.spots <- which(abs(coordinates[,1]-.5) < .15 & abs(coordinates[,2]-.5) < .15)
counts["gene_10",center.spots] <- counts["gene_10",center.spots]+15

x <- haystack_continuous_spatial(coordinates = coordinates, expression = counts,
                                 randomization.count = 50, n.genes.to.randomize = 50)

test_that("haystack spatial works", {
  expect_type(x, "list")
  expect_equal(class(x), "haystack")
  expect_equal(names(x), c("results", "info"))
  expect_equal(class(x$results), "data.frame")
  expect_equal(dim(x$results), c(100, 3))
  expect_equal(x$results["gene_10", "D_KL"], 0.6793742, tolerance = 1e-6)
  expect_equal(x$results["gene_10", "log.p.vals"], -18.20188, tolerance = 1e-6)
})

