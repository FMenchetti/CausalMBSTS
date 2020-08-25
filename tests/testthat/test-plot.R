load("./refmod.Rdata")

test_that("plotting works", {
    x <- sapply( c("impact", "forecast", "ppchecks"),
                function(x) try(plot(causal.1, int.date, type = x, prob = .5)))
    if(interactive()) dev.off()
    expect_true(all(class(x) != 'try-error'))
})
