nbinomTestForMatricesLess <-
function (countsA, countsB, sizeFactorsA, sizeFactorsB, dispsA, 
    dispsB) 
{
    kAs <- rowSums(cbind(countsA))
    kBs <- rowSums(cbind(countsB))
    mus <- rowMeans(cbind(t(t(countsA)/sizeFactorsA), t(t(countsB)/sizeFactorsB)))
    fullVarsA <- pmax(mus * sum(sizeFactorsA) + dispsA * mus^2 * 
        sum(sizeFactorsA^2), mus * sum(sizeFactorsA) * (1 + 1e-08))
    fullVarsB <- pmax(mus * sum(sizeFactorsB) + dispsB * mus^2 * 
        sum(sizeFactorsB^2), mus * sum(sizeFactorsB) * (1 + 1e-08))
    sumDispsA <- (fullVarsA - mus * sum(sizeFactorsA))/(mus * 
        sum(sizeFactorsA))^2
    sumDispsB <- (fullVarsB - mus * sum(sizeFactorsB))/(mus * 
        sum(sizeFactorsB))^2
    sapply(seq(along = kAs), function(i) {
        if (kAs[i] == 0 & kBs[i] == 0) 
            return(NA)
        ks <- 0:(kAs[i] + kBs[i])
        ps <- dnbinom(ks, mu = mus[i] * sum(sizeFactorsA), size = 1/sumDispsA[i]) * 
            dnbinom(kAs[i] + kBs[i] - ks, mu = mus[i] * sum(sizeFactorsB), 
                size = 1/sumDispsB[i])
        pobs <- dnbinom(kAs[i], mu = mus[i] * sum(sizeFactorsA), 
            size = 1/sumDispsA[i]) * dnbinom(kBs[i], mu = mus[i] * 
            sum(sizeFactorsB), size = 1/sumDispsB[i])
        stopifnot(pobs == ps[kAs[i] + 1])
        #if (kAs[i] * sum(sizeFactorsB) < kBs[i] * sum(sizeFactorsA)) 
            numer <- ps[1:(kAs[i] + 1)]
        #else numer <- ps[(kAs[i] + 1):length(ps)]
        min(1,  sum(numer)/sum(ps))
    })
}
