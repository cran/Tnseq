nbinomTestLess <-
function (cds, condA, condB, pvals_only = FALSE, eps = NULL) 
{
    stopifnot(is(cds, "CountDataSet"))
    if (all(is.na(dispTable(cds)))) 
        stop("Call 'estimateDispersions' first.")
    if (dispTable(cds)[condA] == "blind" || dispTable(cds)[condB] == 
        "blind") {
        if (fitInfo(cds, "blind")$sharingMode != "fit-only") 
            warning("You have used 'method=\"blind\"' in estimateDispersion without also setting 'sharingMode=\"fit-only\"'. This will not yield useful results.")
    }
    stopifnot(condA %in% levels(conditions(cds)))
    stopifnot(condB %in% levels(conditions(cds)))
    if (!is.null(eps)) 
        warning("The 'eps' argument is defunct and hence ignored.")
    colA <- conditions(cds) == condA
    colB <- conditions(cds) == condB
    bmv <- getBaseMeansAndVariances(counts(cds)[, colA | colB], 
        sizeFactors(cds)[colA | colB])
    rawScvA <- fData(cds)[, paste("disp", dispTable(cds)[condA], 
        sep = "_")]
    rawScvB <- fData(cds)[, paste("disp", dispTable(cds)[condB], 
        sep = "_")]
    pval <- nbinomTestForMatricesLess(counts(cds)[, colA], counts(cds)[, 
        colB], sizeFactors(cds)[colA], sizeFactors(cds)[colB], 
        rawScvA, rawScvB)
    if (pvals_only) 
        pval
    else {
        bmvA <- getBaseMeansAndVariances(counts(cds)[, colA], 
            sizeFactors(cds)[colA])
        bmvB <- getBaseMeansAndVariances(counts(cds)[, colB], 
            sizeFactors(cds)[colB])
        data.frame(id = rownames(counts(cds)), baseMean = bmv$baseMean, 
            baseMeanA = bmvA$baseMean, baseMeanB = bmvB$baseMean, 
            foldChange = bmvB$baseMean/bmvA$baseMean, log2FoldChange = log2(bmvB$baseMean/bmvA$baseMean), 
            pval = pval, padj = p.adjust(pval, method = "BH"), 
            stringsAsFactors = FALSE)
    }
}
