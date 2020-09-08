context("test-extract_row")

m <- matrix(c(TRUE, TRUE, FALSE, TRUE), nrow = 2) > 0
m <- as(m, "RsparseMatrix")

r <- extract_row_lgRMatrix(m, 1)

test_that("extract_row_lgRMatrix works", {
  expect_equal(r, c(TRUE, FALSE))
})


m <- matrix(1:4, nrow = 2)
m <- as(m, "RsparseMatrix")

r <- extract_row_dgRMatrix(m, 1)

test_that("extract_row_dgRMatrix works", {
  expect_equal(r, c(1, 3))
})
