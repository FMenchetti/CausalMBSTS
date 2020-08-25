test_that("Package examples work", {
    tmp <- tempfile()
    sink(tmp)
    ex <- try(devtools::run_examples('../..'))
    sink()
    unlink(tmp)
    expect_true(class(ex) != 'try-error')
})
