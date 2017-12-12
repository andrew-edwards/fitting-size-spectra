# fitting3repMLEbin-bMinus1.r - testing b=-1 (see Issue #7). 12th December 2017.

# fitting3repMLEbin.r - testing MLEbin, the likelihood method on binned data,
#  to now include in the first manuscript. 3rd May 2016.

# fitting3repBin.r - testing the likelihood method for binned data. So simulating
#  data (as for fitting3rep.r), binning it and then fitting using likelihood
#  based on i) assuming midpoints represent the true values, ii) fitting using
#  likelihood functions that explicitly incorporate the binning.
#  21st October 2015

# fitting3rep.r - as for fitting2rep.r, but now plotting histograms in terms
#  of estimates of b rather than estimates of slope, since they will be clearer.
#  13th July 2015.

rm(list=ls())
print(date())

redo.simulation = TRUE    # Whether or not to redo the simulations
if(!redo.simulation)
  {
  load("fitting3repMLEbin-bMinus1.RData")
  source("../../PLBfunctions.r")
  } else
  {                                                   
  source("../../PLBfunctions.r")
  n = 1000                  # sample size
  b.known = -1              # known fixed value of b
  xmin.known = 1            # known fixed value of xmin
  xmax.known = 1000         # known fixed value of xmax

  num.reps = 10000          # number of times to draw sets of n random numbers.
                            #  (throwing n PLB dice num.reps times)  
  set.seed(42)              # Same seed as for original simulations in
                            #  first manuscript. 
  
  # Record the slope for each method
  
  MLEbin.rep = numeric(num.reps)*NA
         # MLE calculations that explicitly account for binning
  
  # Record the confidence intervals
  MLEbin.rep.conf = data.frame(confMin = MLEbin.rep, confMax = MLEbin.rep)
  
  # Main loop for doing the fitting num.reps times
  for(iii in 1:num.reps)
    {
    if(num.reps > 1000)
      {
      if(iii %in% seq(1000, num.reps, 1000)) print(paste("iii = ", iii))
                                          # show progress
      }
  
    x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)
  
    # Bin the raw data. Linear bins will put almost everything in the first bin
    #  (since a power law), as in Figure 1 of manuscript. So use doubling bins,
    #  which is common for length data in size spectra analyses. 
    
    log2bins.list = log2bins(x)
    num.bins = dim(log2bins.list$binVals)[1]
    # These stopped working in R 3.2.3, need to not have as a table data frame:
    # binBreaks = log2bins.list$binVals[,"binMin"]
    # # Need to append endpoint of final bin:
    # binBreaks = c(binBreaks, as.vector(log2bins.list$binVals[num.bins,
    #    "binMax"]))
    # binCounts = log2bins.list$binVals[,"binCount"]
  
    # This works, but next set of code is shorter:
    #binBreaks = as.data.frame(log2bins.list$binVals[,"binMin"])
    #class(binBreaks) = "data.frame"      # Need to not be a table data frame
    ## Need to append endpoint of final bin:
    #maxOfMaxBin = log2bins.list$binVals[num.bins, "binMax"]
    #class(maxOfMaxBin) = "data.frame"
    #binBreaks = c(binBreaks[,1], maxOfMaxBin[1,1])   # [..] to lose column names
                     # Not the best code, but it works.

    # This is better, and safer since end up as vectors (as original code had):
    binBreaks = log2bins.list$binVals[,"binMin"]$binMin     # Loses column names
    # Need to append endpoint of final bin:
    maxOfMaxBin = log2bins.list$binVals[num.bins, "binMax"]$binMax
    binBreaks = c(binBreaks, maxOfMaxBin)
    
    binCounts = log2bins.list$binVals[,"binCount"]$binCount
    binMids = log2bins.list$binVals[,"binMid"]$binMid  # Midpoints of bins,
                     #  can be used when not properly accounting for the binning
  
    if(!is.integer(binCounts))
      { stop("Need to adapt code for noninteger counts")
      }
                              # will need to do at some point.
     
    xmin = min(binBreaks)
    xmax = max(binBreaks)
  
    
    # Notation:
    # hAAA - h(istrogram) for method AAA.
    
    # MLEbin (maximum likelihood on binned data) calculations.
    # 
    # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
    #  as a starting point for nlm for MLE of b for PLB model.
    # PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1
    
    # Analytical calculation for unbounded power law, to use as a starting point.
    #  If we only have the binned values then need to calculate an approximate:
    sum.log.xBinned2 = sum(log(binMids) * binCounts)
    
    PL.bMLE.binned = 1/( log(min(binBreaks)) -
        sum.log.xBinned2/sum(binCounts)) - 1
  
    # Calculate the MLE for the PLB model using the binned data.
    PLB.bin.minLL = nlm(negLL.PLB.binned, p = PL.bMLE.binned, 
      w = binBreaks, d = binCounts, J = length(binCounts))
  
    PLB.bin.bMLE = PLB.bin.minLL$estimate
    
    MLEbin.rep[iii] = PLB.bin.bMLE
    
    # 95% confidence intervals for MLEbin method.
    PLB.bin.minNegLL = PLB.bin.minLL$minimum
    
    # Values of b to test to obtain confidence interval. For the movement data
    #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
    #  symmetric interval here.
    bvec = seq(PLB.bin.bMLE - 0.5, PLB.bin.bMLE + 0.5, 0.001) 
        
    PLB.bin.LLvals = vector(length=length(bvec))  # -ve log-likelihood for bvec
    for(i in 1:length(bvec))
      {
        PLB.bin.LLvals[i] = negLL.PLB.binned(bvec[i], 
          w = binBreaks, d = binCounts, J = length(binCounts))
      }
    critVal = PLB.bin.minNegLL  + qchisq(0.95,1)/2
                    # 1 degree of freedom, Hilborn and Mangel (1997) p162.
    bIn95 = bvec[ PLB.bin.LLvals < critVal ]
                    # b values in 95% confidence interval
    PLB.bin.MLE.bConf = c(min(bIn95), max(bIn95))
    if(PLB.bin.MLE.bConf[1] == min(bvec) | PLB.bin.MLE.bConf[2] == max(bvec))
      { windows()
      plot(bvec, PLB.bin.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
      }
    MLEbin.rep.conf[iii,] = c(PLB.bin.MLE.bConf[1], PLB.bin.MLE.bConf[2])
  
    }  # End for(iii in 1:num.reps) loop
  
  } # End of if(!redo.simulation) {load("fitting1rep.RData")} else { 

# Prints Latex code for table that summarises the results
# Quartiles:
print("Method & Slope is estimating: & 1st quartile & Median & Mean &
  3rd quartile & Percentage below true \\")
print(paste("MLEbin & $b$ & ", qqtab(MLEbin.rep), "\\hline"))

# 5% and 95% values:
print("Method & Slope represents & 5% quantile & Median & Mean & 95% quantile
  & Percentage below true \\")
print(paste("MLEbin & $b$ & ", qqtab(MLEbin.rep, quants=c(0.05, 0.95)), "\\"))

# Doing all figures in fitting3repMLEbinConf.r

save.image(file = "fitting3repMLEbin-bMinus1.RData")

# To check that final (10,000th) sequence of x values is the same as for
#  fittingrep.RData:

# load("fitting3repMLEbin.RData")
# x[1:10]
# x10bin = x[1:10]
# binRseed = .Random.seed
# load("fitting3rep.RData")
# notbinRseed = .Random.seed
# x10notBin = x[1:10]
# x10bin
# x10notBin
# x10notBin - x10bin
# notbinRseed - binRseed
# All gives 0's, so x values are identical, and would have been throughout
#  the simulations, as expected.
