load("./refmod.Rdata")

test_that("summary works", {
    s <- summary(causal.1)
    expect_true(all(names(s[[1]]) == c("mean", "lower", "upper", "cum.sum", "cum.lower", "cum.upper", "bayes.pval", "pct.causal.eff")))
})
