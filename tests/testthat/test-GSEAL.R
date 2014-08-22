### setup ###

library(org.Sc.sgd.db)

# shared matrix for general use
m = GOMembershipMatrix(org.Sc.sgdGO)

consistency.check = function(m) {
    expect_equal(rowSums(m@matrix != 0), m@geneData$Count)
    expect_equal(colSums(m@matrix != 0), m@colData$Count)
}

context("GO membership matrices")

test_that("We can create a gene set membership matrix", {
    mm = GOMembershipMatrix(org.Sc.sgdGO)
    consistency.check(mm)
    # check that it contains only 0s and ones
    expect_equal(unique(c(as.matrix(mm))), c(0, 1))
})

test_that("min.size and max.size arguments of GOMembershipMatrix work", {
    mm = GOMembershipMatrix(org.Sc.sgdGO, min.size=5, max.size=200)
    consistency.check(mm)
    expect_more_than(min(colSums(mm@matrix)), 4)
    expect_less_than(max(rowSums(mm@matrix)), 201)
})

context("Subsetting a matrix")

test_that("Subsetting a gene matrix leads to correct dimensions", {
    consistency.check(m[2:4, ])
    consistency.check(m[, 50:55])
    consistency.check(m[20:25, 40:50])
    
    expect_equal(dim(m[2:4, ]), c(3, ncol(m)))
    expect_equal(dim(m[, 6:10]), c(nrow(m), 5))
    expect_equal(dim(m[11:50, 20:21]), c(40, 2))

    expect_equal(nrow(m[, 5:7]@colData), 3)
    expect_equal(nrow(m[, 5:7]@geneData), nrow(m@geneData))
    expect_equal(nrow(m[16:20, ]@geneData), 5)
    expect_equal(nrow(m[16:20, ]@colData), nrow(m@colData))
})

test_that("Subsetting a matrix with a numeric vector gets correct rows", {
    # todo
})

test_that("Subsetting a matrix with a logical vector gets correct rows", {
    # todo
})

context("Association testing")

test_that("Wilcoxon rank sum association testing calculates correct p-values", {
    # create a random trait to test
    trait = rnorm(nrow(m))

    # test in all three alternative hypotheses
    testers = sample(which(colSums(m@matrix) > 10), 6)
    for (alternative in c("two.sided", "less", "greater")) {
        m = TestAssociation(m, rownames(m), trait, method="wilcoxon",
                            alternative=alternative)
        manual.tests = apply(as.matrix(m@matrix)[, testers], 2, function(col) {
            wilcox.test(trait[col != 0], trait[col == 0], alternative=alternative)$p.value
        })
        expect_equal(m@colData$pvalue[testers], unname(manual.tests))
    }
})

test_that("T-test association testing calculates correct p-values", {
    # todo: check that it is close enough, may not be exact

    # create a random trait to test
    trait = rnorm(nrow(m))
    
    testers = sample(which(colSums(m@matrix) > 10), 3)
    # test in all three alternative hypotheses
    for (alternative in c("two.sided", "less", "greater")) {
        m = TestAssociation(m, rownames(m), trait, method="t.test",
                            alternative=alternative)
        manual.tests = apply(as.matrix(m@matrix)[, testers], 2, function(col) {
            t.test(trait[col != 0], trait[col == 0], alternative=alternative)$p.value
        })
        expect_equal(m@colData$pvalue[testers], unname(manual.tests))
    }
})

test_that("Hypergeometric association testing calculates correct p-values", {
    # todo
})
