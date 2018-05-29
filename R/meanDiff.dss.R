#' Function to calculate simple mean differences of methylation levels across
#' regions. Takes as input for dmrs a GenomicRanges object, but doesn't 
#' assume that index ranges is present.
meanDiff.dss <- function(bs, dmrs, testCovariate) {
  # convert covariates to column numbers if characters
  if (is.character(testCovariate)) {
    testCovariate <- which(colnames(pData(bs)) == testCovariate)
    if (length(testCovariate) == 0) {
      stop("testCovariate not found in pData(). ", 
           "Please specify a valid testCovariate")
    }
  }
  
  coeff <- seq(2,(2 + length(testCovariate) - 1))
  testCov <- pData(bs)[, testCovariate]
  if (length(unique(testCov)) == 1) {
    message("Warning: only one unique value of the specified ", 
            "covariate of interest.  Assuming null comparison and ", 
            "splitting sample group into two equal groups")
    testCov <- rep(1, length(testCov))
    testCov[seq_len(round(length(testCov)/2))] <- 0
  }
  
  design <- model.matrix(~testCov)
  colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
  
  if (length(unique(design[, coeff])) != 2) {
    message("Not a two-group comparison. Can't compute simple mean ",
            "methylation differences. ", 
            "Returning beta estimates instead")
    return(dmrs$beta)
  } else {
    prop.mat <- getCoverage(bs, type = "M") / 
      getCoverage(bs, type = "Cov")
    levs <- unique(design[, coeff])
    
    ol <- as.data.frame(findOverlaps(dmrs, bs)) %>%
      group_by(queryHits) %>% 
      summarize(s = min(subjectHits),
                e = max(subjectHits))
    miss <- which(!((1:length(dmrs)) %in% ol$queryHits))
    
    indexRanges <- IRanges(ol$s, ol$e)
    prop.mat.dmr <- extractROWS(prop.mat, indexRanges)
    prop.mat1.means <- rowMeans2(prop.mat.dmr[,design[, coeff] == levs[1]],
                                 na.rm=TRUE)
    prop.mat2.means <- rowMeans2(prop.mat.dmr[,design[, coeff] == levs[2]],
                                 na.rm=TRUE)
    meanDiff <- IRanges::mean(IRanges::relist(prop.mat2.means - prop.mat1.means,
                                              indexRanges), na.rm=TRUE)
    md <- rep(NA, length(dmrs))
    if (length(miss)>0){
      md[-miss] <- meanDiff
      return(md)
    }else{
      return(meanDiff)
    }
  }
}
