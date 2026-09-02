test_that("BayesMissingModel rejects a groups vector that does not match values", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)

  expect_error(
    BayesMissingModel(toy$values, toy$groups[-1], cm, parallel = FALSE),
    "one entry per column"
  )
})

test_that("BayesMissingModel rejects fewer than two populated groups", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)
  one_group <- factor(rep("A", ncol(toy$values)), levels = c("A", "B"))

  expect_error(
    BayesMissingModel(toy$values, one_group, cm, parallel = FALSE),
    "at least two non-empty groups"
  )
})

test_that("BayesMissingModel rejects a contrast naming an absent group", {
  toy <- make_toy(c(A = 4, B = 4))
  three <- factor(c(as.character(toy$groups), "C"), levels = c("A", "B", "C"))
  cm <- limma::makeContrasts(contrasts = "C - A", levels = levels(three))

  expect_error(
    BayesMissingModel(toy$values, toy$groups, cm, parallel = FALSE),
    "absent from 'groups'"
  )
})

test_that("BayesMissingModel rejects NA group assignments", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)
  with_na <- toy$groups
  with_na[1] <- NA

  expect_error(
    BayesMissingModel(toy$values, with_na, cm, parallel = FALSE),
    "must not contain NA"
  )
})

test_that("valid unbalanced input passes validation and returns groups as a factor", {
  toy <- make_toy(c(A = 5, B = 3))
  cm <- toy_contrast(toy$groups)

  expect_identical(
    missBayes:::validateGroups(toy$values, as.character(toy$groups), cm),
    factor(as.character(toy$groups))
  )
})
