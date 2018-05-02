# Function to read in methylation count files from GEO and format into 
# BSseq objects. Modified from bsseq function read.bismark.

read.lister <- function (files, sampleNames, rmZeroCov = FALSE, 
                         strandCollapse = TRUE, mc.cores = 1, verbose = TRUE, 
                         BACKEND = NULL, chipBS = FALSE)
{
  if (!is.null(BACKEND) && mc.cores > 1) {
    stop("Currently, 'mc.cores' must be 1 if 'BACKEND' is not NULL")
  }
  if (anyDuplicated(files)) {
    stop("duplicate entries in 'files'")
  }
  if (length(sampleNames) != length(files) || anyDuplicated(sampleNames)) {
    stop("argument 'sampleNames' has to have the same length as argument 'files', without duplicate entries")
  }
  sampleNames <- unname(sampleNames)
  idxes <- seq_along(files)
  names(idxes) <- sampleNames
  allOut <- mclapply(idxes, function(ii) {
    if (verbose) {
      cat(sprintf("[read.bismark] Reading file '%s' ... ",
                  files[ii]))
    }
    ptime1 <- proc.time()
    out <- read.CytosineReportRaw(thisfile = files[ii],
                                  thisSampleName = sampleNames[ii], 
                                  rmZeroCov = rmZeroCov,
                                  BACKEND = BACKEND,
                                  chipBS = chipBS)
    if (strandCollapse) {
      out <- strandCollapse(out)
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
      cat(sprintf("done in %.1f secs\n", stime))
    }
    out
  }, mc.cores = mc.cores)
  if (verbose) {
    cat(sprintf("[read.bismark] Joining samples ... "))
  }
  ptime1 <- proc.time()
  allOut <- combineList(allOut)
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3L]
  if (verbose) {
    cat(sprintf("done in %.1f secs\n", stime))
  }
  allOut
}

read.CytosineReportRaw <- function (thisfile, thisSampleName, rmZeroCov, 
                                    BACKEND = NULL, chipBS = FALSE)
{
  if (isGzipped(thisfile)) {
    stop("File needs to be unzipped and filtered for CGs before proceeding")
  }
  
  out <- fread(thisfile)
  
  if (chipBS){
    if (ncol(out) != 7L) {
      stop("unknown file format")
    }
    gr <- GRanges(seqnames = out[[1]], 
                  ranges = IRanges(start = out[[2]],width = 1), 
                  strand = out[[4]])
    M <- as.matrix(out[[6L]])
    Cov <- as.matrix(out[[7L]])
  }else{
    if (ncol(out) != 6L) {
      stop("unknown file format")
    }
    gr <- GRanges(seqnames = out[[1]], 
                ranges = IRanges(start = out[[2]],width = 1), 
                strand = out[[3]])
    M <- as.matrix(out[[5L]])
    Cov <- as.matrix(out[[6L]])
  }
  
  M <- realize(M, BACKEND = BACKEND)
  Cov <- realize(Cov, BACKEND = BACKEND)
  BSseq(gr = gr, sampleNames = thisSampleName, M = M, Cov = Cov,
        rmZeroCov = rmZeroCov)
}
