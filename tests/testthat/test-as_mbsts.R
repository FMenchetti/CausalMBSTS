load("./refmod.Rdata")

test_that("as.mbsts works", {
    set.seed(12) 
    mbsts.2 <- as.mbsts(y = exd[, c('y1', 'y2', 'y3')], components = c("trend", "seasonal"), seas.period = 7,
                        s0.r = diag(3), s0.eps = diag(3), niter = 10, burn = 1)
    expect_equal(mbsts.1, # from refmod.Rdata
                 mbsts.2)
})
