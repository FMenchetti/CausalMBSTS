load("./refmod.Rdata")

test_that("prediction works", {
    mbsts_pred <- try(predict(mbsts.1, steps.ahead = 10))
    expect_true(class(mbsts_pred) != 'try-error')
})
