test_that("a balanced, group-contiguous run is unchanged by the unbalanced-design fix", {
  # balanced-reference.rds was produced by missBayes 0.99.0, before groups were
  # derived from the group factor rather than from column position. A balanced
  # design whose columns are already grouped together is the regime in which the
  # old positional code was correct, so the fix must be an exact no-op there.
  reference <- readRDS(test_path("balanced-reference.rds"))

  result <- toy_run(make_toy(c(A = 4, B = 4)))

  expect_equal(result, reference, tolerance = 0)
})
