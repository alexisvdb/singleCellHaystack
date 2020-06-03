context("test-get_d_kl")

y <- dat.tsne
param <- singleCellHaystack:::get_parameters_haystack(y[, 1], y[, 2])
ref <- singleCellHaystack:::get_reference(param) # reference distribution
w <- dat.expression > 1 # detection

z <- sapply(1:5, function(k) {
  singleCellHaystack:::get_D_KL(w[k, ], param, ref$Q, ref$pseudo)
})

test_that("multiplication works", {
  expect_equal(z, c(0.06756395, 0.10501541, 0.12736556, 0.12583667, 0.18401621), tolerance = 1e-6)
})
