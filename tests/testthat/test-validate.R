test_that("validate_saint_input accepts a valid list", {
  si <- list(
    inter = data.frame(a = 1),
    prey = data.frame(a = 1),
    bait = data.frame(a = 1)
  )
  expect_true(validate_saint_input(si))
})

test_that("validate_saint_input rejects wrong names", {
  expect_error(validate_saint_input(list(inter = data.frame(), prey = data.frame())),
               "named list")
  expect_error(validate_saint_input(list(a = 1, b = 2, c = 3)),
               "named list")
})

test_that("validate_saint_input requires data frames", {
  si <- list(inter = data.frame(a = 1), prey = data.frame(a = 1), bait = "not a df")
  expect_error(validate_saint_input(si), "data.frame")
})
