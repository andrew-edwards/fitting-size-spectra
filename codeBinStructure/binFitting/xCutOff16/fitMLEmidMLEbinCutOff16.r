# fitMLEmidMLEbinCutOff16.r - updating this (just pathnames) from the 2016 version,
#  rather than the new fitting-size-spectra/ fitMLEbin one, since the 2016
#  version had quite a few changes. 18th May 2018.

# fitMLEmidMLEbinCutOff16.r - testing MLEmid and MLEbin methods with xmin=1 and
#  a cut-off of 16. Need to change sample size also, and check how that gets
#  used. 8th August 2016.

# fitMLEmidMLEbin.r - testing MLEmid and MLEbin methods with xmin=16.
#  20th July 2016.

# fitMLEmidMLEbin.r - testing MLEmid and MLEbin methods. 3rd June 2016.

# fitting3repMLEbin.r - testing MLEbin, the likelihood method on binned data,
#  to now include in the first manuscript. 3rd May 2016.

rm(list=ls())
print(date())
start = proc.time()
redo.simulation = TRUE    # Whether or not to redo the simulations
if(!redo.simulation)
  {
  load("fitMLEmidMLEbinCutOff16.RData")
  source("../../../code/PLBfunctions.r")
  source("../../countsFunctions.r")     # repeated here for developing
  } else
  {
  source("../../../code/PLBfunctions.r")
  source("../../countsFunctions.r")

  n.desired = 1000          # sample size desired after cut off.
  b.known = -2              # known fixed value of b
  xmin.known = 1            # known fixed value of xmin. Make it a power of 2
  if(!is.wholenumber(log2(xmin.known)) | !is.wholenumber(xmin.known))
       { stop("If want xmin.known to not be an integer and not be an integer
                 power of 2 then need to edit binData; may just need to think
                 about whether can just set startInteger = FALSE, but will
                 need to edit binData to be able to define the binWidth for
                 2k method. Will need to think about this for real data; maybe
                 best to just remove the first bin (and a fit a range that is
                 encompassed by the binned data).")
      }
  xmax.known = 1000         # known fixed value of xmax

  cutOff = 16               # cut off: not using data < cutOff.
                            # **IF CHANGE CUT OFF and re-run then fix x> issue.
  full.sample.size = 1.5 * n.desired /
                (1 - pPLB(cutOff, b = b.known, xmin=xmin.known, xmax=xmax.known))
                            # full.sample.size to draw from to expect to get n
                            # values >cutOff.
                            # Originally did exactly this and then used the
                            #  values >= cutOff, which gave a distribution of
                            #  realised sample sizes, minimum of which was 858
                            #  (x.min=1, cutOff=16, b=-2). So 1000/858 = 1.165
                            #  suggests multipling full.sample.size by 1.2 should
                            #  get 1000 for each sample, so use 1.5 just to be
                            #  safe (this isn't the time consuming part of the
                            #  code). This number will change for different
                            #  parameter values so may need to be calculated
                            #  analytically.

  num.reps = 10000          # number of times to draw sets of n random numbers.
                            #  (throwing n PLB dice num.reps times)
  set.seed(42)              # Same seed as for original simulations in
                            #  first manuscript.

# From Rowan, to create empty 3-d array with names. I wouldn't name the rows.


  binType = list(1, 5, 10, "2k")   # or use substring("lin1", 4) or grep
                                   #  and then adapt binData.
                                   # Must be numeric, for linear bins, or "2k"
                                   #  so define as a list
  binTypes = length(binType)
  binType.name = binType
  binType.name[!binType %in% "2k"] = paste0("Linear ",
                  binType[which(!binType %in% "2k")])

  MLEmethod.name = c("MLEmid", "MLEbin")    # If these change then need to change likelihood
  MLEmethods = length(MLEmethod.name)       #  calls below, and table output; not automatic.

  # Do an array for MLE's and then an array for confMin and confMax
  MLE.array = array(NA, dim=c(num.reps, binTypes, MLEmethods),
      dimnames=list(1:num.reps, unlist(binType.name), MLEmethod.name))
       # no need to name the rows, just index by simulation number
  # MLE.array[i,j,k] is random sample i, bin type j, MLE method j

  # Record the confidence intervals
    MLEconf.array = array(NA, dim=c(num.reps, binTypes, MLEmethods, 2),
      dimnames=list(1:num.reps, unlist(binType.name), MLEmethod.name,
      c("confMin", "confMax")))
  # MLEconf.array[i,j,k, ] is confidence interval [c(confMin, confMax)] for
  #   random sample i, bin type j, MLE method k

  #  Record the realised sample size
  #  n.realised = vector(length=num.reps)

  # Main loop for doing the fitting num.reps times
  for(i in 1:num.reps)
  {
    if(num.reps > 1000)
      {
      if(i %in% seq(1000, num.reps, 1000)) print(paste("i = ", i))
                                          # show progress
      }

    x.full = rPLB(full.sample.size, b = b.known, xmin = xmin.known,
        xmax = xmax.known)
    x.cut = x.full[x.full >= cutOff]
    if(length(x.cut) < n.desired)
        {
           save.image(file = "fitMLEmidMLEbinCutOff16temp.RData")
           stop("Need to increase full.sample.size; try and continue
                   simulations if possible.")
        }
    x = x.cut[1:n.desired]
    #n.realised[i] = length(x)    # Don't think n is actually needed

    for(j in 1:binTypes)                    # Loop over binning type
    {
      bins.list = binData(x, binWidth=binType[[j]])  # 1, 2, 5 or "2k", etc.
      num.bins = dim(bins.list$binVals)[1]

      binBreaks = bins.list$binVals[,"binMin"]$binMin   # Loses the column names
      maxOfMaxBin = bins.list$binVals[num.bins, "binMax"]$binMax
      binBreaks = c(binBreaks, maxOfMaxBin) # Append endpoint of final bin

      binCounts = bins.list$binVals[,"binCount"]$binCount
      binMids = bins.list$binVals[,"binMid"]$binMid  # Midpoints of bins

      if(sum(!is.wholenumber(binCounts)) > 0)
        { stop("Need to adapt code for noninteger counts")
        }
                                # will need to do at some point. Though I think
                                #  negLL.PLB.counts can handle it, but don't
                                #  worry about until everything more
                                #  functionalised.

      sumCntLogMids = sum(binCounts * log(binMids))

      # MLEmid (maximum likelihood using midpoints) calculations.

      # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
      #  as a starting point for nlm for MLE of b for PLB model.
      PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts)) - 1

      # Note that using the min and max of binBreaks for xmin and xmax,
      #  because this method acknowledges that data are binned (but
      #  negLL.PLB.counts can act on just counts of discrete values). The
      #  min and max of binBreaks are the MLE's of xmin and xmax.
      MLEmid.res = calcLike(negLL.PLB.counts, p=PL.bMLE, x=binMids,
        c=binCounts, K=num.bins, xmin=min(binBreaks), xmax=max(binBreaks),
        sumclogx=sumCntLogMids)

      MLE.array[i,j,"MLEmid"] = MLEmid.res$MLE
      MLEconf.array[i,j,"MLEmid", ] = MLEmid.res$conf

      # MLEbin (maximum likelihood on binned data) calculations.

      MLEbin.res = calcLike(negLL.PLB.binned, p=PL.bMLE,
          w = binBreaks, d = binCounts, J = length(binCounts))

      MLE.array[i,j,"MLEbin"] = MLEbin.res$MLE
      MLEconf.array[i,j,"MLEbin", ] = MLEbin.res$conf

    }   # End of for(j in 1:binTypes) loop

  }  # End of for(i in 1:num.reps) loop

} # End of if(!redo.simulation) {load("fitting1rep.RData")} else {

# Prints Latex code for table that summarises the results


# Doing all figures in fitting3repMLEbinConf.r

write.csv(bins.list$binVals, "binValsCutOff16.csv")  # row.names=FALSE - put in at some point
write.csv(bins.list$indiv, "indivCutOff16.csv")
save.image(file = "fitMLEmidMLEbinCutOff16.RData")

print(date())
print("Time to run program, seconds:")
print(proc.time()-start)

