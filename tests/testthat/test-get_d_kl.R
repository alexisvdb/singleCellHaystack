context("test-get_d_kl")

set.seed(123)
x <- scale(mtcars)
y <- Rtsne::Rtsne(x, perplexity = 5)$Y
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])
ref <- singleCellHaystack:::get_reference(param) # reference distribution
w <- x > 1 # detection

z <- sapply(1:5, function(k) {
  singleCellHaystack:::get_D_KL(w[k, ], param, ref$Q, ref$pseudo)
})

test_that("multiplication works", {
  expect_equal(z, c(0.5109014, 0.5109014, 0.4967910, 0.4584567, 0.3187935), tolerance = 1e-6)
})
