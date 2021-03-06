context("Wrapper functions")

test_that("zeroes do not interfere with execution", {
  expect_message(svsample(c(0, y), draws = draws, burnin = burnin, quiet = TRUE), "Argument 'y' *")
})

test_that("svsample executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

test_that("svsample with nu executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

test_that("svsample with rho executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

