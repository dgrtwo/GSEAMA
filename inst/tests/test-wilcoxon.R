test_that("Wilcoxon rank-sum test generates correct p-values", {
    data(factorial)
    mm = MembershipMatrix(organism = "Sc.sgd", ontology = "BP", min.size = 5, max.size = 250)
    
    # perform three hypothesis tests, one for each alternative
    for (ah in c("two.sided", "greater", "less")) {
        wilcoxon.mm = TestEnrichment(mm, factorial$ORF, factorial$RNA.Seq.logFC, method = "wilcoxon", alternative=ah)
        
        # check that it is internally consistent
        expect_that(NROW(wilcoxon.mm@setData), equals(NCOL(wilcoxon.mm@matrix)))
        expect_that(NROW(wilcoxon.mm@geneData), equals(NROW(wilcoxon.mm@matrix)))
        
        # check the first 6 and the last 6 pvalues
        check.vals = as.numeric(c(1:6, tail(seq_along(wilcoxon.mm@setData$pvalue))))
        for (ind in check.vals) {
            b = (wilcoxon.mm@matrix[, ind] > 0)
            manual.pval = wilcox.test(wilcoxon.mm@geneData$y[b], wilcoxon.mm@geneData$y[!b], alternative=ah)$p.value
            expect_that(wilcoxon.mm@setData$pvalue[ind], equals(manual.pval))
            # check that hypothesis test means what we want it to
            r = rank(wilcoxon.mm@geneData$y)
            members.greater = mean(r[b]) > mean(r[!b])
            if (ah == "greater") {
                expect_that(wilcoxon.mm@setData$pvalue[ind] < .5,
                            equals(members.greater))
            }
            if (ah == "less") {
                expect_that(wilcoxon.mm@setData$pvalue[ind] < .5,
                            equals(!members.greater))
            }
        }    
    }
})

