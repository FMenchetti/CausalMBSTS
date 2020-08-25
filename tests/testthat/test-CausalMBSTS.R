
load("./refmod.Rdata")

test_that("CausalMBSTS works", {
    set.seed(12)
    causal.2 <- CausalMBSTS(exd[, c('y1', 'y2', 'y3')], components = c("trend", "seasonal"),
                            seas.period = 7, X = exd[, c('x1', 'x2', 'x3', 'x4')], dates = dates,
                            int.date = int.date, s0.r = 0.1*diag(3), s0.eps = 0.1*diag(3),
                            niter = 10, burn = 1, horizon = as.Date(c('2019-02-19','2019-03-10')))
    expect_equal(causal.1, ## from refmod.Rdata
                 causal.2)
})

test_that("CausalMBSTS works with data.frame input", {
    exd <- as.data.frame(exd)
    causal.3 <- CausalMBSTS(exd[, c('y1', 'y2', 'y3')], components = c("trend", "seasonal"),
                            seas.period = 7, X = exd[, c('x1', 'x2', 'x3', 'x4')], dates = dates,
                            int.date = int.date, s0.r = 0.1*diag(3), s0.eps = 0.1*diag(3),
                            niter = 10, burn = 1, horizon = as.Date(c('2019-02-19','2019-03-10')))
    expect_equal(causal.1, ## from refmod.Rdata
                 causal.3)
})
