# countsFunctions.r - functions for dealing with binned data and associated
#  counts in each bin. 7th June 2016.

require(plotrix)       # for axis.break function
require(dplyr)

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
    {
    # Is x a whole number? Taken from ?is.integer and added documentation.
    # Args:
    #  x: vector
    #  tol: tolerance
    # Returns:
    #  vector of TRUE/FALSE corresponding to each element of x
    abs(x - round(x)) < tol
    }


negLL.PLB.counts = function(b, x, c, K=length(c), xmin=min(x), xmax=max(x),
    sumclogx = sum(c * log(x)))
  {
  # Calculates the negative log-likelihood of the parameters b, xmin and xmax
  #  given count data for the PLB model. Returns the negative log-likelihood.
  #  Will be called by nlm or similar, but xmin and xmax will be estimated
  #  as the min of lowest bin and max of the largest, no need for any numerics.
  #  For testing the MLEmid methods (using midpoints of bins), then
  #  give xmin and xmax explicitly as the lowest and highest bin breaks
  #  because the x values correspond to bins. But if x just represents
  #  counts of discrete values then no need to specify xmin and xmax, they
  #  will be automatically determined as min(x) and max(x), respectively,
  #  although it can be good to specify them to avoid repeated calculation.    
  #
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood.
  #   x: vector of length K corresponding to data values x_k, with a
  #    corresponding count c being the number of times that xk is repeated.  
  #   c: vector of length K giving the counts c_k for each k=1, 2, 3, ..., K.
  #       Must have c[1]>0 and c[K]>0, i.e. non-zero counts for first and last
  #       x_k. Note that the c_k do not have to be integer-valued.    
  #   K: number of c_k (length of c).
  #   xmin: minimum value of x_k, as an input to avoid repeatedly calculating.
  #   xmax: maximum value of x_k, as an input to avoid repeatedly calculating.
  #   sumclogx: sum( c * log(x) ), to avoid repeatedly calculating.
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
    if(xmin <= 0 | xmin >= xmax | length(x) != K | length(c) != K |
         c[1] == 0 | c[K] == 0 | min(c) < 0)
         stop("Parameters out of bounds in negLL.PLB.counts")
    n = sum(c)
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumclogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumclogx
      }
    return(neglogLL)
  }


LBNbiom.method.counts = function(valCounts, binBreaks = NULL, lowerCutOff = 0)
    {
    # Use the log-binning with normalisation technique to calculate
    #  the slope of the biomass size spectra, for count data.
    #  Slope is from fitting
    #  a linear regression of log10(normalised biomass in bin)
    #  against log10(midpoint of bin). Bins can be defined by user,
    #  else are created to double in size. Also calculates slope
    #  for biomasses not being normalised. 
    # 
    # Args:
    #  valCounts: data.frame (and can be tbl_df) with columns bodyMass
    #   and Number (which is the count for each body mass). bodyMass can
    #   represent midpoints, say, of existing bins, or be the actual
    #   species-specific converted-to-bodyMass. Number can be non-integer,
    #   which can arise from standardising, say, trawl data to be per hour.
    #  binBreaks: breaks for the bins to be used to bin the data and
    #   then fit the regression. If not provided then it calculates
    #   them as bin widths that double in size that encompass the data,
    #   resulting in binBreaks   ..., 0.25, 0.5, 1, 2, 4, 8, 16,....
    #   as necessary.
    #  lowerCutOff: body mass value representing the lower cut off for
    #   the range being fit.
    #
    # Returns:
    #  list containing:
    #   valCounts2: dataframe valCounts with extra columns binMin, the
    #    minimum of the bin into which that bodyMass falls, and biomass,
    #    the biomass corresponding to bodyMass * Number.
    #   binVals: dataframe with a row for each bin, where the
    #    columns are: binMid, binMin, binMax, binWidth - midpoint, minimum,
    #                  maximum, and width, respectively, of the bin
    #                 totalBiom - total biomass in that bin
    #                 totalBiomNorm - normalised total biomass in that bin,
    #                  defined as totalBiom / binWidth
    #                 log10.... - log10 of some of the above quantities
    #  norm.lm: lm() result of the linear regression fit using normalised
    #    biomass in each bin
    #  norm.slope: slope of the linear regression fit using normalised
    #    biomass in each bin
    #  unNorm.lm: lm() result of the linear regression fit when not
    #    normalising the biomass in each bin
    #  unNorm.slope: slope of the linear regression fit when not
    #    normalising the biomass in each bin
    # 
        require(dplyr)
        if(!is.data.frame(valCounts))
            { stop("valCounts not a data.frame in LBNbiom.method.counts")}
        if(anyNA(valCounts))
            { stop("valCounts contains NA's in LBNbiom.method.counts") }
        if(min(valCounts$bodyMass) <= 0)
            { stop("valCountsbodyMass needs to be >0 in LBNbiom.method.counts") }
        #   x = bodyMass
        xmin = min(valCounts$bodyMass)
        xmax = max(valCounts$bodyMass)
        #
        if(is.null(binBreaks))
           {
            binBreaks = 2^( floor(log2(xmin)) : ceiling(log2(xmax)) )
           } else
           { 
            if(min(diff(binBreaks)) < 0) { stop("binBreaks need to be increasing
                                              in LBNbiom.method.counts")}
           }
        if(!(lowerCutOff %in% c(0, binBreaks))) { stop("need lowerCutOff
                           to be 0 or one of the binBreaks in LBNbiom.method.counts") }
           # We could allow a user to specify a lowerCutoff that was not
           #  a binBreak, but it just makes defining the binBreaks a bit more
           #  fiddly -- code could be modified if a user wished. Although in
           #  practice people would plot the binned data and then choose which
           #  points (binned counts) to ignore when fitting the regression.
        valCounts2 = mutate(valCounts,
            binMin = binBreaks[findInterval(bodyMass, binBreaks,
                rightmost.close=TRUE)], biomass = bodyMass * Number)
                                        # need total biomass for each row
        if(max(valCounts2$binMin)  == binBreaks[length(binBreaks)])
           { stop("check binBreaks in LBNbiom.method.counts") }   # Shouldn't occur

        binVals = summarise(group_by(valCounts2, binMin), binCount = sum(Number),
            totalBiom = sum(biomass))
        # No guarantee that every binMin value hLlBNbiom.mids shows up here
        missing1 = setdiff(binBreaks[-length(binBreaks)], binVals$binMin)
        missing = cbind(binMin = missing1, binCount = rep(0, length(missing1)),
            totalBiom = rep(0, length(missing1)))
        binVals = rbind(binVals, missing)
        binVals = tbl_df(binVals)
        binVals = arrange(binVals, binMin)
        if( max( abs( binBreaks[-length(binBreaks)] - binVals$binMin) ) > 0)
            { stop("check binVals in LBNbiom.method.counts") }
        # So all the breaks except final show up as binMin, and have been ordered.
        #  Therefore to add binMax just do:
        binVals = mutate(binVals, binMax = binBreaks[-1], binWidth = binMax - binMin,
            binMid = binMin + binWidth/2,  totalBiomNorm = totalBiom / binWidth )
        binVals = mutate(binVals, log10binMid = log10(binMid),
            log10totalBiom = log10(totalBiom),
            log10totalBiomNorm = log10(totalBiomNorm))
        binVals[ is.infinite(binVals$log10totalBiom),
                  "log10totalBiom"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[ is.infinite(binVals$log10totalBiomNorm),
                  "log10totalBiomNorm"] = NA
        binVals = mutate(binVals, aboveCutOff = (binMid > lowerCutOff))
                  # aboveCutOff is TRUE/FALSE for the regression
        unNorm.lm = lm( log10totalBiom ~ log10binMid,
            data = filter(binVals, aboveCutOff),
            na.action=na.omit)
        unNorm.slope = unNorm.lm$coeff[2]
        unNorm.conf = confint(unNorm.lm, "log10binMid", 0.95)
        
        norm.lm = lm( log10totalBiomNorm ~ log10binMid,
            data = filter(binVals, aboveCutOff),
            na.action=na.omit)
        norm.slope = norm.lm$coeff[2]
        norm.conf = confint(norm.lm, "log10binMid", 0.95)        
        y = list(valCounts2 = valCounts2, binVals = binVals,
            unNorm.lm = unNorm.lm, unNorm.slope = unNorm.slope,
            unNorm.conf = unNorm.conf,
            norm.lm = norm.lm, norm.slope = norm.slope, norm.conf = norm.conf,
            lowerCutOff = lowerCutOff)
        return(y)
       }

Llin.method.counts = function(valCounts, num.bins = NULL, binBreaks = NULL)
    {
    # The Llin method, which is plotting binned counts on log-linear
    #  axes and then fitting a regression, for count data.
    # 
    # Args:
    #  valCounts: data.frame (and can be tbl_df) with columns bodyMass
    #   and Number (which is the count for each body mass). bodyMass can
    #   represent midpoints, say, of existing bins, or be the actual
    #   species-specific converted-to-bodyMass. Number can be non-integer,
    #   which can arise from standardising, say, trawl data to be per hour.
    #  num.bins: number of bins to be used, though this is only a suggestion
    #   since can get over-ridden by hist(). Need to specify num.bins OR
    #   binBreaks.
    #  binBreaks: breaks for the bins to be used to bin the data and
    #   then fit the regression.
    #
    # Returns:
    #  list containing:
    #   mids: midpoint of bins
    #   log.counts: log(counts) in each bin
    #   counts: counts in each bin
    #   lm: results of the linear regression of log.counts ~ mids
    #   slope: slope of the linear regression fit
    #   breaks: bin breaks
    #   confVals: 95% confidence interval of the fitted slope
    #   stdError: standard error of the fitted slope
        require(dplyr)   
        if(!is.data.frame(valCounts))
            { stop("valCounts not a data.frame in Llin.method.counts")}
        if(anyNA(valCounts))
            { stop("valCounts contains NA's in Llin.method.counts") }
        if(min(valCounts$bodyMass) <= 0)
            { stop("valCountsbodyMass needs to be >0 in Llin.method.counts") }
        if(is.null(binBreaks) & is.null(num.bins))
            { stop("need binBreaks or num.bins in Llin.method.counts") }
        if(!is.null(binBreaks) & !is.null(num.bins))
            { stop("need one of binBreaks OR num.bins in Llin.method.counts") }
        if(!is.null(binBreaks))      # use to get the bin breaks
           { 
            if(min(diff(binBreaks)) < 0) stop("binBreaks need to be increasing")
            hLlin.temp = hist(valCounts$bodyMass, breaks = binBreaks, plot=FALSE)
                      # will give error if binBreaks don't span bodyMass values
           }  else    # { breaks = num.bins }
           {
            hLlin.temp = hist(valCounts$bodyMass, breaks = num.bins, plot=FALSE)
           }          # Still lets hist select the breaks; num.bins is a guide
        breaks = hLlin.temp$breaks
        hLlin.mids = hLlin.temp$mids               

        valCounts2 = mutate(valCounts,
            binMid = hLlin.mids[findInterval(bodyMass, breaks,
                rightmost.close=TRUE)])
        binVals = summarise(group_by(valCounts2, binMid), binCount = sum(Number))
        # No guarantee that every binMid hLlin.mids shows up here (unlikely,
        #  given linear binning and long-tailed data). Even though need to
        #  exclude 0 counts (since get logged) from lm, still best to return
        #  values for all bins, especially since the bins can be specified.
        missing1 = setdiff(hLlin.mids, binVals$binMid)
        missing = cbind(binMid = missing1, binCount = rep(0, length(missing1)))
        binVals = rbind(binVals, missing)
        binVals = tbl_df(binVals)
        binVals = arrange(binVals, binMid)
        if( max( abs( hLlin.mids - binVals$binMid) ) > 0)
            { stop("check binVals in Llin.method.counts") }
        hLlin.log.counts = log(binVals$binCount)
                  # binVals$binCount was hLlin$counts
        hLlin.log.counts[ is.infinite(hLlin.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

        hLlin.lm = lm( hLlin.log.counts ~ hLlin.mids, na.action=na.omit)

        y = list(mids = hLlin.mids, log.counts = hLlin.log.counts,
            counts = binVals$binCount, lm = hLlin.lm, slope = hLlin.lm$coeff[2],
            breaks = breaks,
            confVals = confint(hLlin.lm, "hLlin.mids",0.95),
            stdErr = coef(summary(hLlin.lm))["hLlin.mids", "Std. Error"])
        return(y)
       }


eightMethods.count = function(data = data, oneYear = 1980,
  figName = "eightMCnt" )
  {
  # Computes exponents for a dataset using all eight methods, explicitly
  #  using the counts (number of each bodyMass), which can be non-integer.
  # Args: 
  #  oneYear: the year of data to use, from that in the multiple years contained
  #   in data.
  #  data: local data frame that has a unique row for every combination
  #   of Year, SpecCode and LngtClass. The `Number' column is
  #   the number of observed individuals of that species in that length
  #   class in that year. `bodyMass' is the body mass representative of such
  #   an individual, as calculated previously by LWa * LngtClass ^ LWb.
  #  figName: figure name, will get appended by -Year for each year, to create
  #   a .png for each year.
  # Returns:
  #   data frame with one row for each method, with columns Year, Method,
  #    b (estimate of b from that method), confMin (lower end of 95% confidence
  #    interval of b for that method), and confMax (upper end of 95% confidence
  #    interval of b for that method).
  #   plots a .png figure paste(figName, "-", oneYear, ".png") of the 
  #    fits for each of the eight methods. Does .png instead of .eps since
  #    .eps was >150Mb for 1980 data, .png is 66Kb.
  #
  # 
  # In eightMethods(), had to expand data to get a vector, x, of individual
  #  fish sizes (lengths or weights), which is how the original methods
  #  functions are written.
  # Now adding explicit methods in eightMethods.counts, such as
  #  Llin.method.counts, to deal explicitly with counts, and that
  #  should also work for non-integer counts. For integer
  #  counts, the expansion to give x should give the same results.
  dataForYear = filter(data, Year == oneYear)
  valCounts = select(dataForYear, bodyMass, Number)  # Just Number of individuals
                            # with each bodyMass (Number can be non-integer).
                            # Don't worry about species here.
  valCounts = ungroup(valCounts)      # Else it retains group info, affecting
                                      #  results for, e.g., findInterval(),
                                      #  so safer to ungroup.
  # x = rep(dataForYear$bodyMass, dataForYear$Number)
              # 3.3 million for 1980 for old nSea15, but it was quick

  # May want some of these:
  # log.x = log(x)                     
  # sum.log.x = sum( log.x ) 
  xmin = min(valCounts$bodyMass)
  xmax = max(valCounts$bodyMass)

  figheight = 7 # 5.6     For 4x2 figure
  figwidth = 5.7    # 5.7 inches for JAE

  num.bins = 8   # number of bins for standard histogram and Llin method, though
                 #  this is only a suggestion (and can get overridden). Daan used
                 #  8 bins.

  # postscript("nSea1980-fitting2a.eps", height = figheight,
  #           width = figwidth, horizontal=FALSE, paper="special")
  # Postscript plots all 3million points, ends up >150Mb file. So try .png:
  png(paste(figName, "Cnt-", oneYear, ".png", sep=""), height = figheight,
             width = figwidth, res=300,, units="in")
  
  par(mfrow=c(4,2))
  oldmai = par("mai")    #  0.95625 0.76875 0.76875 0.39375  inches I think,
                         #   think may be indpt of fig size
  par(mai=c(0.4, 0.5, 0.16, 0.3))  # Affects all figures if don't change again
                                  #  Changed 3rd from 0.05 to 0.13
  mgpVals = c(1.6,0.5,0)            # mgp values   2.0, 0.5, 0
  
  # Notation:
  # hAAA - h(istrogram) for method AAA.
  
  # Llin method - plotting binned data on log-linear axes then fitting regression
  #  as done by Daan et al. 2005.
  # hLlin.list = Llin.method(x, num.bins = num.bins)
  hLlin.list = Llin.method.counts(valCounts, num.bins = num.bins)

  #eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
  #    "b" = hLlin.list$slope,
  #    "confMin"= hLlin.list$confVals[1], "confMax"= hLlin.list$confVals[2]) 
  # That gives hLlin.mids as the name of the row, so does this though
  hLlin.b = hLlin.list$slope
  hLlin.confMin = hLlin.list$confVals[1]
  hLlin.confMax = hLlin.list$confVals[2]
  hLlin.stdErr = hLlin.list$stdErr
  eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
      "b" = hLlin.b, "confMin"= hLlin.confMin, "confMax"= hLlin.confMax,
      "stdErr" = hLlin.stdErr,
      row.names=NULL)

  plot( hLlin.list$mids, hLlin.list$log.counts, 
     xlab=expression(paste("Bin midpoints for data ", italic(x))),
     ylab = "Log (count)", mgp=mgpVals)   # xlim=c(0, 400),

  lm.line(hLlin.list$mids, hLlin.list$lm)
  inset = c(0, -0.04)     # inset distance of legend
  legend("topright", paste("(a) Llin slope=", signif(hLlin.list$slope, 3)),
       bty="n", inset=inset)

  mtext(
   paste("                                                           ",
         oneYear) )
  # LT method - plotting binned data on log-log axes then fitting regression,
  #  as done by Boldt et al. 2005, natural log of counts plotted against natural
  #  log of size-class midpoints.

  # Use Llin method's binning.
  hLT.log.mids = log(hLlin.list$mids)
  hLT.log.counts = log(hLlin.list$counts)
  hLT.log.counts[ is.infinite(hLT.log.counts) ] = NA  
                    # lm can't cope with -Inf, which appear if 0 counts in a bin
  
  hLT.lm = lm( hLT.log.counts ~ hLT.log.mids, na.action=na.omit)
  hLT.slope = hLT.lm$coeff[2]
  hLT.conf = confint(hLT.lm, "hLT.log.mids", 0.95)
  hLT.stdErr = coef(summary(hLT.lm))["hLT.log.mids", "Std. Error"]
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LT"),
      "b" = hLT.slope, "confMin"= hLT.conf[1], "confMax"= hLT.conf[2], 
      "stdErr" = hLT.stdErr, row.names=NULL))
   
  plot( hLT.log.mids, hLT.log.counts, 
       xlab=expression(paste("Log (bin midpoints for data ", italic(x), ")")), 
       ylab = "Log (count)", mgp=mgpVals)

  lm.line(hLT.log.mids, hLT.lm)
  legend("topright", paste("(b) LT b=", signif(hLT.slope, 3)), bty="n",
         inset=inset)

  # LTplus1 method - plotting linearly binned data on log-log axes then fitting
  #  regression of log10(counts+1) vs log10(midpoint of bins), as done by
  #  Dulvy et al. (2004).

  # Use Llin method's binning.
  hLTplus1.log10.mids = log10(hLlin.list$mids)
  hLTplus1.log10.counts = log10(hLlin.list$counts + 1)
  hLTplus1.log10.counts[ is.infinite(hLTplus1.log10.counts) ] = NA  
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
                  #  but the + 1 avoids this issue here
  hLTplus1.lm = lm( hLTplus1.log10.counts ~ hLTplus1.log10.mids, 
      na.action=na.omit)
  hLTplus1.slope = hLTplus1.lm$coeff[2]
  hLTplus1.stdErr = coef(summary(hLTplus1.lm))[
      "hLTplus1.log10.mids", "Std. Error"]  

  hLTplus1.conf = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95)
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LTplus1"),
      "b" = hLTplus1.slope, "confMin"= hLTplus1.conf[1], 
      "confMax"= hLTplus1.conf[2], "stdErr" = hLTplus1.stdErr, row.names=NULL))
  
  plot( hLTplus1.log10.mids, hLTplus1.log10.counts, 
       xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")), 
       ylab = "Log10 (count+1)", mgp=mgpVals)
 
  lm.line(hLTplus1.log10.mids, hLTplus1.lm)
  legend("topright", paste("(c) LTplus1 b=", signif(hLTplus1.slope, 3)),
       bty="n", inset=inset)

  # LBmiz method - binning data using log10 bins, plotting results on natural
  #  log axes (as in mizer). Mizer does abundance size spectrum or biomass
  #  size spectrum - the latter multiplies abundance by the min of each bin
  #  (see below).

  # Construction of bins is as follows, from Finlay
  #  Scott:
  # The bins dimensions can be specified by the user by passing in min_w, max_w
  #  [min values for the lowest and highest bins] and no_w arguments [number of
  #  bins]. These are then used:
  #    w <- 10^(seq(from=log10(min_w), to=log10(max_w), length.out=no_w))
  #    dw <- diff(w)
  #    dw[no_w] <- dw[no_w-1] # Set final dw as same as penultimate bin
  #
  #The w values are the break points of the bins (the start of the bin).
  # Regarding the regression, x and w will have the same length since x
  # is just the abundance (numbers or biomass) at size w.
  # From mizer's getCommunitySlopeCode.r:
  #  "Calculates the slope of the community abundance through time by
  #  performing a linear regression on the logged total numerical abundance
  #  at weight and logged weights (natural logs, not log to base 10, 
  #  are used)."  So regress log(total counts) against log(weights)
  #  (not log10 and not normalised). And it's actually
  #   on the minima of the bins (their w).
 
  hLBmiz.num.bins = num.bins
  
  beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est
  
  # hLBmiz.bins = c(beta^(0:(k-1)) * xmin, xmax)
  hLBmiz.bins = c(beta^(0:(hLBmiz.num.bins-1)) * xmin, xmax)
  hLBmiz.mins = hLBmiz.bins[-length(hLBmiz.bins)]    # min of each bin
   # Mizer bin specification, with final bin being same width as penultimate bin

# Adapting from Llin.method.counts:
  LBmiz.valCounts = mutate(valCounts,
                binMin = hLBmiz.mins[findInterval(bodyMass, hLBmiz.bins,
                rightmost.closed=TRUE)])    # Note that this would retain groups
  LBmiz.binVals = summarise(group_by(LBmiz.valCounts, binMin),
                binCount = sum(Number))
        # No guarantee that every binMid hLlin.mids shows up here (unlikely,
        #  given linear binning and long-tailed data). Even though need to
        #  exclude 0 counts (since get logged) from lm, still best to return
        #  values for all bins, especially since the bins can be specified.
  LBmiz.missing1 = setdiff(hLBmiz.mins, LBmiz.binVals$binMin)
  LBmiz.missing = cbind(binMin = LBmiz.missing1,
      binCount = rep(0, length(LBmiz.missing1)))
  LBmiz.binVals = rbind(LBmiz.binVals, LBmiz.missing)
                    # works even if LBmiz.missing has no rows
  LBmiz.binVals = tbl_df(LBmiz.binVals)
  LBmiz.binVals = arrange(LBmiz.binVals, binMin)
  if( max( abs( hLBmiz.mins - LBmiz.binVals$binMin) ) > 0)
            { stop("check LBmiz.binVals$binMin in eightMethods.counts") }


#  hLlin.log.counts = log(binVals$binCount)
                  # binVals$binCount was hLlin$counts
#        hLlin.log.counts[ is.infinite(hLlin.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.log.min.of.bins = log(LBmiz.binVals$binMin)    # min of bins
  hLBmiz.log.counts = log(LBmiz.binVals$binCount)
  hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA  
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)

  # hLBmiz.slope = hLBmiz.lm$coeff[2]
  # Need to subtract 1, since want to work in terms of b not slopes now
  hLBmiz.b = hLBmiz.lm$coeff[2] - 1

  hLBmiz.conf = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95) - 1
  hLBmiz.stdErr = coef(summary(hLBmiz.lm))[
      "hLBmiz.log.min.of.bins", "Std. Error"]
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBmiz"),
      "b" = hLBmiz.b, "confMin"= hLBmiz.conf[1], 
      "confMax"= hLBmiz.conf[2], "stdErr" = hLBmiz.stdErr, row.names=NULL))

  plot( hLBmiz.log.min.of.bins, hLBmiz.log.counts,
     xlab=expression(paste("Log (minima of bins for data ", italic(x), ")")),
     ylab = "Log (count)", mgp=mgpVals)
     # axes=FALSE, xaxs="i", yaxs="i", , xlim=c(log10(1), log10(650)), 
     #  ylim=c(log10(0.7), log10(1100)))  # So axes are logged

  lm.line(hLBmiz.log.min.of.bins, hLBmiz.lm)
  legend("bottomleft", paste("(d) LBmiz b=", signif(hLBmiz.b, 3)), 
         bty="n", inset=inset)

  # LBbiom method - binning data using log2 bins, calculating biomass not counts
  #  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
  #  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

  hLBNbiom.list = LBNbiom.method.counts(valCounts)    # Does this method and the next.

  hLBbiom.b = hLBNbiom.list[["unNorm.slope"]] - 2
  hLBbiom.conf = hLBNbiom.list[["unNorm.conf"]] - 2
  hLBbiom.stdErr = coef(summary(hLBNbiom.list$unNorm.lm))[
      "log10binMid", "Std. Error"]
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBbiom"),
      "b" = hLBbiom.b, "confMin"= hLBbiom.conf[1], 
      "confMax"= hLBbiom.conf[2], "stdErr" = hLBbiom.stdErr, row.names=NULL))
  
  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiom,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["unNorm.lm"]])
  legend("bottomleft", paste("(e) LBbiom b=",
     signif(hLBbiom.b, 3)), bty="n", inset=c(-0.08, -0.04))

  # LBNbiom method - on biomass, not counts, as per Julia's 2005 paper.
  #  log2 bins of bodymass, sum the total biomass in each bin, normalise
  #  biomasses by binwidths, fit regression to log10(normalised biomass) v
  #  log10(midpoint of bin).
  
  # hLBNbiom.list = LBNbiom.method(x) - already done above

  hLBNbiom.b = hLBNbiom.list[["norm.slope"]] - 1
  hLBNbiom.conf = hLBNbiom.list[["norm.conf"]] - 1
  hLBNbiom.stdErr = coef(summary(hLBNbiom.list$norm.lm))[
      "log10binMid", "Std. Error"]
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBNbiom"),
      "b" = hLBNbiom.b, "confMin"= hLBNbiom.conf[1], 
      "confMax"= hLBNbiom.conf[2], "stdErr" = hLBNbiom.stdErr, row.names=NULL))
  
  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])
  legend("bottomleft", paste("(f) LBNbiom b=",
     signif(hLBNbiom.b, 3)), bty="n", inset=inset)

  # Cumulative Distribution, LCD method
  # logProp = log((1:length(x))/length(x))               # x equivalent:
  # This method should really split up the cumProp into, say, 1000 values
  #  since for the regression each point gets weighted the same but this
  #  isn't quite right. But this is how someone would likely do it.
  #  To plot the results for the MLE method below I'm generating 1000 values.
  LCD.valCounts = select(valCounts, bodyMass, Number)
  LCD.valCounts = arrange(LCD.valCounts, desc(bodyMass))
                                        # x.sorted = sort(x, decreasing=TRUE)
  sumNumber = sum(LCD.valCounts$Number)
  LCD.valCounts = mutate(LCD.valCounts, cumSum = cumsum(Number))
                                        # 1:length(x)
  LCD.valCounts = mutate(LCD.valCounts, cumProp = cumSum / sumNumber)
                                        # logProp = log((1:length(x))/length(x))

  LCD.valCounts = mutate(LCD.valCounts, logBodyMass = log(bodyMass),
                            logCumProp = log(cumProp))
                                       # logSorted = log(x.sorted)
  hLCD.lm = lm(logCumProp ~ logBodyMass, data = LCD.valCounts)
                                        # hLCD.lm = lm(logProp ~ logSorted)
  hLCD.slope = hLCD.lm$coeff[2]

  hLCD.b = hLCD.lm$coeff[2] - 1
  hLCD.conf = confint(hLCD.lm, "logBodyMass", 0.95) - 1
  hLCD.stdErr = coef(summary(hLCD.lm))[
      "logBodyMass", "Std. Error"]
 
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LCD"),
      "b" = hLCD.b, "confMin"= hLCD.conf[1], 
      "confMax"= hLCD.conf[2], "stdErr" = hLCD.stdErr, row.names=NULL))

   plot(LCD.valCounts$logBodyMass, LCD.valCounts$logCumProp, main="",
     xlab=expression(paste("Log ", italic(x))),
     ylab=expression( paste("Log (prop. of ", values >= italic(x), ")")),
     mgp=mgpVals) # , axes=FALSE)
     #xlim=c(0.8, 1000), xaxs="i", ylim=c(0.0008, 1), yaxs="i",

  lm.line(LCD.valCounts$logBodyMass, hLCD.lm, col="red")
  # murankfreq = 1 - fitsortedlog10$coeff[2]       # mu = 1 - slope
  legend("bottomleft", paste("(g) LCD b=", signif(hLCD.b, 3)), bty="n",
       inset=inset)

  
  # MLE (maximum likelihood method) calculations.
  MLE.valCounts = select(valCounts, bodyMass, Number)
  # Should be faster to group repeated values, just in case that's not done:  
  MLE.valCounts = summarise(group_by(valCounts, bodyMass), Count = sum(Number))
  MLE.valCounts = arrange(MLE.valCounts, desc(bodyMass))
                                        # x.sorted = sort(x, decreasing=TRUE)
  sumCounts = sum(MLE.valCounts$Count)
  if(abs( sumCounts - sumNumber) > 0.001) 
      { stop("Check sumCounts in eightMethods.count()") }
  MLE.K = dim(MLE.valCounts)[1]         # Number of bodyMass values
  if(MLE.valCounts[1, "Count"] == 0 |
     MLE.valCounts[MLE.K, "Count"] == 0)
      { stop("Need first and last counts to be zero in
          eightMethods.count() for MLE method")}
  # Can adapt the code, but better to just to take such zero counts out of data
  # Adapting (for counts) analytical value of MLE b for PL model
  #  (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  MLE.xmin = min(MLE.valCounts$bodyMass)
  MLE.xmax = max(MLE.valCounts$bodyMass)
  MLE.sumCntLogMass = sum(MLE.valCounts$Count * log(MLE.valCounts$bodyMass)) 
  PL.bMLE.counts = 1/( log(MLE.xmin) - MLE.sumCntLogMass/sumCounts) - 1
      
  PLB.minLL =  nlm(negLL.PLB.counts, p=PL.bMLE.counts,
      x=MLE.valCounts$bodyMass, c=MLE.valCounts$Count, K=MLE.K,
      xmin=MLE.xmin, xmax=MLE.xmax, sumclogx=MLE.sumCntLogMass)
                                        #, print.level=2 )
  
  PLB.bMLE = PLB.minLL$estimate
  
  # 95% confidence intervals for MLE method.
  
  PLB.minNegLL = PLB.minLL$minimum
  
  # Values of b to test to obtain confidence interval. For the real movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.
  
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.00001) 

      
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
          PLB.LLvals[i] = negLL.PLB.counts(bvec[i], x=MLE.valCounts$bodyMass,
           c=MLE.valCounts$Count, K=MLE.K,
           xmin=MLE.xmin, xmax=MLE.xmax, sumclogx=MLE.sumCntLogMass)
      }
  critVal = PLB.minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
                      # b values in 95% confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLB.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
    }

  #  Assume for now that confidence interval is PLB.bMLE +/- 1.96 * stdErr
  #   just to quickly calc stdErr. So stdErr = -(confMin - PLB.bMLE)/1.96.
  #      PLB.bMLE - PLB.MLE.bConf - is symmetric anyway, implying quadratic
  #      likelihood profile, implying normal approximation is okay. So use
  #      it now to go backwards, to properly calculate should do the
  #      Fisher information. stdErr = 1 / (sqrt( d^2 logLik /db^2 at MLE))
  #      though that is fiddly to derive.   
  PLB.MLE.stdErr = mean(abs((PLB.MLE.bConf - PLB.bMLE)/1.96))
  
  # MLE.rep.xmax[iii] = xmax   - not storing xmax for now
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("MLE"),
      "b" = PLB.bMLE, "confMin"= PLB.MLE.bConf[1], 
      "confMax"= PLB.MLE.bConf[2], "stdErr" = PLB.MLE.stdErr, row.names=NULL))
  # Can delete: TO HERE  ***INteresting - how to plot with the repeated values. If integers, then plot will look different depending if use counts or not. For integers should expand out. Maybe the way to go, and for LCD also, is for non-integer counts to create 1,000 individual body masses values that 'represent' the observed distribution. Fiddly but doable, but will be the best representation, and will show up repeated values......HERE***
  # To plot rank/frequency style plot, because of non-integer counts want to
  #  calculate cumulative proportions (as in LCD), and then generate 1,000
  #  individual body masses values that 'represent' the observed distribution,
  #  and then plot y axis as Proportion >= x, with 1,000 points plotted.
  MLE.valCounts = mutate(MLE.valCounts, cumSum = cumsum(Count))
                                        # x equivalent:  1:length(x)
  MLE.valCounts = mutate(MLE.valCounts, cumProp = cumSum / sumCounts)
                                        # logProp = log((1:length(x))/length(x))
  # So if you had a sample of 1,000 that includes xmin and xmax, what would
  #  the remaining 998 body masses be?:
  if(MLE.valCounts[dim(MLE.valCounts)[1],"cumProp"] != 1)
      { stop("Check MLE.valcounts in eightMethods.counts()")  }
  MLE.sim = tbl_df(data.frame(cumPropSim =
      seq(MLE.valCounts[1,"cumProp"], 1, length=ceiling(sumCounts))))  
                                        # simulated cumulative proportions
                                        # or maybe do on sumSum, though props
                                        #  better capture endpoints I think?
                                        # length does ceiling anyway, so make
                                        #  explicit here. Did try just 1000
                                        #  example fish, but that doesn't capture
                                        #  the details in the tail (since
                                        #  actual sample size is 33,592.02), so
                                        #  this does what a sample of 33,593
                                        #  would likely look like.
  MLE.sim = mutate(MLE.sim, bodyMassSim = MLE.valCounts[
                findInterval(cumPropSim, MLE.valCounts$cumProp), "bodyMass"])

      # findInterval():
      #  vec = MLE.valCounts$cumProp     # breakpoints
      # Find interval containing each of the  elements of
      #  x = MLE.sim$cumPropSim  (x in findInterval() terminology)
      
#        LCD.valCounts = mutate(LCD.valCounts, logBodyMass = log(bodyMass),
      #                      logCumProp = log(cumProp))
                                       # logSorted = log(x.sorted)
  plot(MLE.sim$bodyMassSim, MLE.sim$cumPropSim, log="xy",
       xlab=expression(paste("Values, ", italic(x))),
       ylab=expression( paste("Proportion of ", values >= x)), mgp=mgpVals)

  x.PLB = seq(xmin, xmax, length=1000)     # x values to plot PLB
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                    xmax = max(x.PLB))) # * sumCounts
  lines(x.PLB, y.PLB, col="red") #, lty=5)
  
  legend("bottomleft", paste("(h) MLE b=", signif(PLB.bMLE, 3)), bty="n",
         inset=inset)
  
  # To add the curves at the limits of the 95% confidence interval:
  #for(i in c(1, length(bIn95)))   # for(i in 1:length(bIn95))  to see all vals
  #    {
  #      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
  #                  xmax = max(x.PLB))) * length(x), col="red", lty=3)
  #    }  
  
  dev.off()
  return(eightMethodsRes)
}

# Function to plot and fit time series of exponent as estimated by
#  a given method. Called from nSea15analysis.Snw
timeSerPlot = function(bForYears, legName, method, weightReg = FALSE,
    bCol="black",
    pchVal = 20, cexVal = 1, confCol="black", confThick = 1,
    xLim = NULL, yLim = NULL, xLab="",
    yLab = expression(paste("Estimate of ", italic(b)), sep=""),
    xTicksSmallInc = NULL, xTicksSmallTck = 0.01,
    yLabels = TRUE, yTicksSmallInc = NULL, yTicksSmallTck = 0.01,
    legPos = "topleft", newPlot = TRUE,
    regPlot = TRUE,
    regColNotSig = "darkgrey", regColSig = "red",
    legExtra = NULL, legExtraPos = "topleft", legExtraCol = "",
    insetVal = c(-0.08, -0.06))         # insetVal2 = c(-0.08, 0.07)) 
                                        # yTicks = seq(0, 300, 50), 
    {
    # Plotting time series of estimated b with confidence intervals, for
    #  data that are analysed year-by-year by a single method. And then
    #  fit a linear regression with ** intervals.
    #  This will get called eight times to produce a comparison figure of
    #  the methods. 
    # Args:
    #  bForYears: data-frame with columns Year, Method (optional), b
    #   (the estimate of b), confMin and confMax (the 95% lower and upper
    #   confidence limits) and stdError (the standard error of the estimate of b).
    #  legName: legend name for that panel
    #  method: method used to obtain the inputted estimates of b
    #  weightReg: TRUE if doing weighted regression (using standard errors)
    #   or FALSE to not do weighted regression.
    #  bCol: colour for points for b
    #  pchVal: pch for points for b
    #  cexVal: size of points for b
    #  confCol: colour for confidence intervals for b
    #  confThick: thickness of vertical line for confidence intervals
    #  xLim: xlim range
    #  yLim: limits for y-axis
    #  yLab: label for y axis
    #  xTicksSmallInc: increments for where to have small (unlabelled)
    #   tickmarks on x-axis
    #  xTicksSmallTck: tick length for small (unlabelled) tickmarks on x-axis
    #  yTicksSmallInc: increments for where to have small (unlabelled)
    #   tickmarks on y-axis
    #  yTicksSmallTck: tick length for small (unlabelled) tickmarks on y-axis
    #  yLabels: whether or not to label main tickmarks on y-axis
    #  legPos: legend position
    #  newPlot: TRUE to create a new plot, FALSE to add to existing
    #  regPlot: TRUE to plot the regression line and conf intervals    
    #  regColNotSig: colour for regression line (and its confidence intervals)
    #   if the trend is not significant
    #  regColSig: colour for regression line (and its confidence intervals)
    #   if the trend is significant
    #  legExtra: extra manually-specified legend (e.g. to distinguish two
    #   sets of results)
    #  legExtraPos: position for extra manually-specified legend
    #  legExtraCol: colours (vector) for extra manually-specified legend        
    #  insetVal: inset shift for naming the panel
    #  # insetVal2: inset shift for printing observed coverage percentage
    if(is.null(xLim))
        {
            #rangeVal = range(repConf)
            #xLim = c(floor(rangeVal[1]), ceiling(rangeVal[2]))
            # min(repConf[,1])), ceiling(max(repConf[,2])))
            xLim = range(bForYears$Year)
        }
    if(is.null(yLim))        # just do the yLim for this set of results
        {
        #yLim = c(floor(min(bForYears$confMin, na.rm=TRUE)), 
        #  ceiling(max(bForYears$confMax, na.rm=TRUE)))
            # **Need pretty here?:
        yLim = range(c(bForYears$confMin, bForYears$confMax), na.rm=TRUE)
        }        # For Llin and LT can get only two bins and so NaN for
                 #  conf intervals, so need na.rm.
    if(newPlot)
       {
       plot(bForYears$Year, bForYears$b, xlim=xLim, ylim=yLim, col=bCol,
         pch=pchVal, cex=cexVal, xlab=xLab, ylab=yLab)    # yaxt="n")
    
       legend(legPos, legName, bty="n", inset=insetVal)
       if(!is.null(yTicksSmallInc))
           { yTicksSmall = seq(yLim[1], yLim[-1], by=yTicksSmallInc)
             axis(2, at = yTicksSmall, labels = rep("", length(yTicksSmall)),
                tck=-yTicksSmallTck)
           }  
       if(!is.null(xTicksSmallInc))
           { xTicksSmall = seq(xLim[1], xLim[-1], by=xTicksSmallInc)
             axis(1, at = xTicksSmall, labels = rep("", length(xTicksSmall)),
                tck=-xTicksSmallTck)
           }  
       # Confidence intervals (instead of plotCI from regress2.Snw):
       segments(x0=bForYears$Year, y0=bForYears$confMin, x1=bForYears$Year,
             y1=bForYears$confMax)
       if(!is.null(legExtra)) legend(legExtraPos, legExtra, bty="n",
                                     col=legExtraCol, pch=pchVal, cex=cexVal)
       } else    # Add to existing plot
       {
       points(bForYears$Year, bForYears$b, col=bCol,
         pch=pchVal, cex=cexVal)
       segments(x0=bForYears$Year, y0=bForYears$confMin, x1=bForYears$Year,
             y1=bForYears$confMax)
       }
 
           
    # Now just fit a linear regression through the points,
    #  and colour code it red if significant trend and grey if not. This
    #  is taken and adapted from regress2.Snw from RBR14 assessment.
    if(weightReg == TRUE)
        { lm = lm(b ~ Year, data = bForYears, weights = 1/(stdErr^2))   } else
        { lm = lm(b ~ Year, data = bForYears) }
          
    yearInc = seq(xLim[1], xLim[2], 0.1)
    p.conf = predict(lm, newdata=data.frame(Year=yearInc), interval="confidence")
    pVal = summary(lm)$coeff["Year",4]
    if(regPlot)
      {
        if (pVal > 0.05) regCol = regColNotSig else regCol= regColSig
        lm.line(xLim, lm, col=regCol)
        matlines(yearInc, p.conf[ ,c("lwr", "upr")], col=regCol, 
                   lty=2)
      }
    confVals = confint(lm, "Year", level=0.95)   

    res = data.frame(Method = method,
       Low = confVals[1], Trend = lm$coeff[2], 
       High = confVals[2], p = pVal, Rsquared = summary(lm)$r.squared,
       adjRsquared = summary(lm)$adj.r.squared, row.names=NULL) 
    return(res)
}

# Copying log2bins from PLBfunctions.r to here, renaming it
#  binData, then generalising it to also use linear bins of any size.
#  It is called from fitMLEmidMLEbin.r for now, so
binData = function(x = NULL, counts = NULL, binWidth = NULL, binBreaks = NULL,
                   startInteger = TRUE)
    {
    #  Construct bins that start from floor(min(x)) or min(x) and either double
    #  in size or are of equal width, and encompass the data.
    # 
    # Args:
    #  x: vector of individual values (e.g. body masses).
    #   OR
    #  counts: dataframe (or array) with first column being an x value
    #  (e.g. body mass), and second column being the counts of the
    #   number of individuals for that value.
    #   Only x or counts can be specified.
    # binWidth: type of bins to use:
    #      "2k" will result in binBreaks that with:
    #         - startInteger=TRUE are powers of 2, i.e.
    #            ..., 0.25, 0.5, 1, 2, 4, 8, 16,.... 
    #         - with startInteger=FALSE are bins that double in size and
    #            start with min(x); not yet implemented, since have to
    #            think about what the width of the first bin should be.    
    #      value 'a' will result in binBreaks are separated by a and span the
    #       data, that with:
    #         - startInteger=TRUE start from z = floor(min(x)) and are then
    #             z, z+a, z+2a, z+3a, ....   (if z = 0 then power-law cannot
    #             be so need to use startInteger=FALSE)
    #         - startInteger=FALSE start from z = min(x) and are then
    #             z, z+a, z+2a, z+3a, ....   
    # binBreaks: pre-defined bin breaks as a vector. Only binWidth
    #   or binBreaks can be specified.
    # startInteger: TRUE or FALSE, whether to start the binBreaks at an integer
    #   power of 2 (for method "2k") or an integer. See binWidth above,
    #   startInteger is ignored if binBreaks is specified.
    # Returns:
    #  list containing:
    #   
    #   indiv: dataframe with a row for each x value, where the
    #    columns are: x - original x values.
    #                 binMid, binMin, binMax, binWidth - midpoint, minimum,
    #                  maximum, and width, respectively, of the bin within
    #                  which the x value falls.
    #     If indiv has >=10^6 rows then it isn't saved.
    #     If counts was specified then an equivalent x
    #      vector is created and is column x (i.e. x values are repeated). May
    #      not be the most efficient way, but it easiest to program.
    #     May not need indiv returned, but it needs to be calculated anyway.
    #
    #   binVals: dataframe with a row for each new bin, where the
    #    columns are:
    #      binMid, binMin, binMax, binWidth - midpoint, minimum,
    #         maximum, and width, respectively, of the bin
    #      binCount - total number of individuals in that bin
    #      binCountNorm - binCount / binWidth
    #      binSum - sum of individual values in that bin
    #           (appropriate if x represents biomass, but not length)
    #      binSumNorm - binSum / binWidth    
    #      log10.... - log10 of some of the above quantities
        require(dplyr)
        if(!is.null(x) & !is.null(counts)) {
            stop("need only one of x or counts in binData") }
        if(is.null(x) & is.null(counts)) {
            stop("need x or counts in binData") }
        if(!is.null(x)) {
          if(!is.vector(x))stop("x not a vector in binData")
          if(anyNA(x)) stop("x contains NA's in binData")
          if(min(x) <= 0)stop("x needs to be >0 in binData")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in binData")
          if(min(counts[,1]) < 0) {
              stop("x values in counts need to be >= 0 in binData") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in binData") }
          if(sum(!is.wholenumber(counts[,2])) > 0) {
              stop("numbers in counts need to be integers in binData;
                      for non-integer count see a new function. Currently,
                      a new function has no name. Or it may be easier to
                      adapt binData.") }
          } 
        if(is.null(binWidth) & is.null(binBreaks)) {
            stop("need one of binWidth or binBreaks in binData") }
        if(!is.null(binWidth) & !is.null(binBreaks)) {
            stop("need only one of binWidth or binBreaks in binData") }
        if(startInteger != "TRUE" & startInteger != "FALSE"){
            stop("startInteger must be TRUE or FALSE in binData") }
        # As for LBNbiom.method(), could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code. But do this
        #  at some point for noninteger counts.
        if(is.null(x))  
           { x = rep(counts[,1], counts[,2])
             minx = min(counts[,1])
             maxx = max(counts[,1])
           } else
           { minx = min(x)
             maxx = max(x)
           }  
        #
        if(!is.null(binBreaks))
           {
           if(minx < min(binBreaks) | maxx > max(binBreaks) )
             { stop("binBreaks do not span data in binData")
             }           
           } else           # create binBreaks based on binWidth
           {
           if(binWidth == "2k")
             {
             if(startInteger)
               { binBreaks = 2^( floor(log2(minx)) : ceiling(log2(maxx)) )
               } else
               { stop("startInteger currently needs to be TRUE when
                   binWidth = 2k")
               }  
             } else     # If not "2k"
             {
             if(!is.numeric(binWidth))
               { stop("binWidth must be 2k or a number (in quotes is okay
                         in quotes) in binData().")
               }
             # startInteger says whether to start from an integer value or
             #  start from min(x),
             z = floor(minx) * startInteger + minx * !startInteger
             binBreaks = seq( z, by=binWidth,
                        length = ceiling( (maxx - z)/binWidth) + 1)
             }
           }

        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))        
        indiv$binMin =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-1])
        #
        indiv$binMid = as.numeric(as.character(indiv$binMid))
        indiv$binMin = as.numeric(as.character(indiv$binMin))
        indiv$binMax = as.numeric(as.character(indiv$binMax))
        # Now calculate biomass in each bin class:
        binVals = summarise(group_by(indiv, binMid),
            binMin = unique(binMin),
            binMax = unique(binMax),
            binWidth = binMax - binMin,
            binCount = length(x),
            binCountNorm = binCount / binWidth,
            binSum = sum(x),
            binSumNorm = binSum / binWidth )
           # binWidth uses new columns binMax and binMin
        # Indices for minima of bins that have zero counts and so do not
        #  appear in binVals yet:
        emptyBinMinInd = !(binBreaks[-length(binBreaks)] %in% binVals$binMin)
        emptyBinMin = binBreaks[emptyBinMinInd]
        empties = length(emptyBinMin)
        emptyBinMax = binBreaks[-1][emptyBinMinInd]
        emptyBinWidth = emptyBinMax - emptyBinMin
        emptyBinMid = emptyBinMin + emptyBinWidth/2

        emptyVals = as.data.frame(cbind(emptyBinMid, emptyBinMin,
            emptyBinMax, emptyBinWidth, matrix(0, nrow=empties, ncol=4)))
        names(emptyVals) = names(binVals)
        binVals = rbind(binVals, emptyVals)         # still a local df
        
        binVals = binVals[order(binVals$binMid),]   # order by binMid

        binVals = mutate(binVals, log10binMid = log10(binMid),
            log10binCount = log10(binCount),
            log10binSum = log10(binSum),
            log10binSumNorm = log10(binSumNorm))
        binVals[ is.infinite(binVals$log10binCount),
                  "log10binCount"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[ is.infinite(binVals$log10binCountNorm),
                  "log10binCountNorm"] = NA
        binVals[ is.infinite(binVals$log10binSum),
                  "log10binSum"] = NA
        binVals[ is.infinite(binVals$log10binSumNorm),
                  "log10binSumNorm"] = NA
        if(dim(indiv)[1] < 10^6) {       # only save indiv if not too big
          y = list(indiv = indiv, binVals = binVals)
          } else
          {    
          y = list(binVals = binVals)
          }
        return(y)
       }


profLike = function(negLL.fn, MLE, minNegLL, vecDiff=0.5, vecInc=0.001, ...)
    {
    # Profile log likelihood method to calculate 95% confidence interval
    #  for a given negative log-likelihood function, MLE and minimum of
    #  of the negative log-likelihood function.
    # Args:
    #   negLL.fn: negative log-likelihood function that take arguments
    #    (parameters and data) in ... and returns a negative
    #    log-likelihood value.
    #   MLE: MLE, as already calculated.
    #   minNegLL: the minimum of the negative log-likelihood, at the MLE.
    #   vecDiff: the range over which to test the negative log-likelihood
    #    to construct the confidence interval. Default is 0.5 and a symmetric
    #    range is tested for fitting size spectra, since for movement data
    #    sets in Table 2 of Edwards (2011) the intervals were symmetric, so
    #    symmetric seems likely.
    #   vecInc: increments to try, the accuracy of the resulting bounds
    #    will depend on this. Note that a resulting interval of, say,
    #    (-2.123, -1.987) means that that interval is contained within the
    #    'true' 95% interval, which is itself contained within (-2.124, -1.986).
    #    The 'true' bounds lie between the stated lower bounds and between 
    #    the stated upper bounds. So reduce vecInc if further accuracy is needed.
    # Returns:
    #  two-component vector of the 95% confidence interval.
    vec = seq(MLE - vecDiff, MLE + vecDiff, vecInc)
                 # Values of parameter to test to obtain confidence interval

    # LLvals = vector(length=length(bvec))
    LLvals = sapply(X=vec, FUN=negLL.fn, ...)
    critVal = minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
    vecIn95 = vec[ LLvals < critVal ]
                      # values in 95% confidence interval
    conf = c(min(vecIn95), max(vecIn95))
    if(conf[1] == min(vec) | conf[2] == max(vec))
      { windows()
        plot(vec, LLvals)
        abline(h = critVal, col="red")
        stop("Need to make vec larger - see R window")   # Could automate
      }
    return(conf)
}


calcLike = function(negLL.fn, p, vecDiff=0.5, vecInc=0.001, ...)
    {
    # Calculate the maximum likelihood estimate and the 95% confidence
    #  interval using the profile log-likelihood method, for a given
    #  negative log-likelihood function and arguments (parameters and data)
    #  for the negative log-likelihood function.
    # Args:
    #   negLL.fn: negative log-likelihood function that take arguments
    #    (parameters and data) in ... and returns a negative
    #    log-likelihood value.
    #   p: starting point to calculate the maximum likelihood estimate.
    #   vecDiff: the range over which to test the negative log-likelihood
    #    to construct the confidence interval. Default is 0.5 and a symmetric
    #    range is tested for fitting size spectra, since for movement data
    #    sets in Table 2 of Edwards (2011) the intervals were symmetric, so
    #    symmetric seems likely.
    #   vecInc: increments to try, the accuracy of the resulting bounds
    #    will depend on this. Note that a resulting interval of, say,
    #    (-2.123, -1.987) means that that interval is contained within the
    #    'true' 95% interval, which is itself contained within (-2.124, -1.986).
    #    The 'true' bounds lie between the stated lower bounds and between 
    #    the stated upper bounds. So reduce vecInc if further accuracy is needed.
    # Returns:
    #  list containing:
    #    MLE: the maximum likelihood estimate
    #    conf: the 95% confidence interval of the MLE.
    minLL = nlm(f=negLL.fn, p=p, ...)
   
    MLE = minLL$estimate
    conf = profLike(negLL.fn=negLL.fn, MLE=MLE, minNegLL=minLL$minimum, ...)
    res = list(MLE = MLE, conf = conf)
    return(res)
}


totalBiomass = function(bvec, r = NULL, xmin=NULL, xmax=NULL, n=1000)
  {
    # Calculate the totalBiomass (in same units as xmin or nondimensionalised
    #  scaled to xmin) for given values of b (in the vector bvec), n and
    #  either both xmin and xmax or just r.
    # 
    # Args:
    #  bvec: vector of size-spectrum exponents.
    #  r: ratio of xmax/xmin. Need to specificy r or xmin and xmax, all scalars. 
    #  xmin: minimum allowable body size.
    #  xmax: maximum allowable body size.
    #  n: number of individuals.
      
    #
    # Returns:
    #  vector of total biomass values corresponding to the values of b in
    #   bvec; total biomass has same units as xmin (if xmin and xmax specified),
    #   or is nondimensional if r is specified.
    #  
    if( !( !is.null(xmin) & !is.null(xmax) & is.null(r) |
           is.null(xmin) & is.null(xmax) & !is.null(r) )) {
            stop("need either r or xmin and xmax in totalBiomass") }
    if(is.null(r))
      {
        returnDiml = TRUE      # if TRUE then return dimensional T
        r = xmax/xmin
      } else
      {
        returnDiml = FALSE
      }

    Tvec = n * (bvec+1)/(bvec+2) * (r^(bvec+2) - 1) / (r^(bvec+1) - 1)
    Tvec[bvec == -1] = n * (r - 1) / log(r)
    Tvec[bvec == -2] = n * r * log(r) / (r - 1)

    if(returnDiml)
      {
        output = Tvec * xmin
      } else
      {
        output = Tvec
      }
    return(output)
  }  

negLL.PLB.binned.species = function(b, dataBinForLike, dataBinForLikeSummary)
  {
  # Calculates the negative log-likelihood of b for the PLB model,
  #  given binned data where the bins can be different for each species.
  #  Returns the negative log-likelihood. **USES ANDY'S ORIGINAL LIKELIHOOD
  #   FUNCTION -- SEE negLL.PLB.bins.species() to use Mike's simpler one.
  #  Will be called by nlm or similar, but xmin and xmax will just be estimated
  #  as the min of lowest bin and max of the largest bin (that is their MLEs),
  #  no need to do numerically. See Appendix of second manuscript for derivation.
  #
  # ***COPIED FROM negLL.PLB.binned for now, for which b=-1 needs correcting
  #    
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   dataBinForLike: table data frame where each row is the count in a bin
  #      of a species, where the columns [and corresponding mathematical
  #      notation in Appendix] are:
  #          SpecCode: code for each species [s]
  #          wmin: lower bound of the bin [w_{sj}, where j is the bin number]
  #          wmax: upper bound of the bin [w_{s,j+1}]
  #          Number: count in that bin for that species [d_{sj}].
  #      For each species the first and last bins must be non-empty, i.e.
  #      w_{s1}, w_{s,J_s +1} > 0. **Write code to check that before
  #      calling this function (since this gets repeatedly called).
  #   dataBinForLikeSummary: table data frame with one row for each species,
  #      giving
  #      the minimum lower bound [w_{s1}] and maximum upper bound [w_{s,J_s +1}]
  #      and the number of counts for that species, where J_s is the number of
  #      bins for species s (won't need to explicilty specify). Columns are:
  #          SpecCode: code for each species [s]
  #          wminSpecies: minimum lower bound [w_{s1}]
  #          waxSpecies: maximum upper bound [w_{s,J_s +1}]
  #          n_s: total number of counts for species s [n_s]
  #
  #    
  #   negLL.PLB.binned:
  #   w: vector of length J+1 giving the bin breaks w_1, w_2, ..., w_{J+1}
  #   d: vector of length J giving the count in each bin. Must have d_1, d_J >0
  #   J: number of bins (length of d)
  #   xmin: minimum value of bins, as an input to avoid repeatedly calculating
  #   xmax: maximum value of bins, as an input to avoid repeatedly calculating
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
  # **MOVE THESE TO PRE-PROCESS function:
  #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
  #       d[1] == 0 | d[J] == 0 | min(d) < 0)
  #       stop("Parameters out of bounds in negLL.PLB")
  if(b != -1)
      {  # First component of (A.55) and then sum:
         temp1 = mutate(dataBinForLikeSummary, comp1 = - n_s *
             log( abs( wmaxSpecies^(b+1) - wminSpecies^(b+1) ) ) )
         comp1Sum = sum(temp1$comp1)      
         #
         temp2 = mutate(dataBinForLike, comp2 = Number *
             log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
         comp2Sum = sum(temp2$comp2)
         #
         neglogLL = - comp1Sum - comp2Sum     # *Negative* log-likelihood
      } else
      { stop("NOT DONE b=-1 yet; adapt, but first correct, negLL.PLB.binned")
      }
    return(neglogLL)
  }


negLL.PLB.bins.species = function(b, dataBinForLike, n, xmin, xmax)
  {
  # Calculates the negative log-likelihood of b for the PLB model,
  #  given binned data where the bins can be different for each species.
  #  Returns the negative log-likelihood, calculated using Mike's simpler
  #  likelihood function.
  #  Will be called by nlm or similar, but xmin and xmax will just be estimated
  #  as the min of lowest bin and max of the largest bin (that is their MLEs),
  #  no need to do numerically. See Appendix of second manuscript for derivation.
  #
  # ***COPIED FROM negLL.PLB.binned.species, for which b=-1 may need correcting
  #    
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   dataBinForLike: table data frame where each row is the count in a bin
  #      of a species, where the columns [and corresponding mathematical
  #      notation in Appendix] are:
  #          SpecCode: code for each species [s]
  #          wmin: lower bound of the bin [w_{sj}, where j is the bin number]
  #          wmax: upper bound of the bin [w_{s,j+1}]
  #          Number: count in that bin for that species [d_{sj}].
  #   n: total number of counts [n = \Sum_{sj} d_{sj} over all s and j]
  #   xmin: maximum likelihood estimate for xmin [xmin = min_{sj} w_{s,1} ]   
  #   xmax: maximum likelihood estimate for xmax [xmax = max_{sj} w_{s,J_s+1} ]
  #
  #
  #   
  #    THINK NOT NEEDED in negLL.PLB.bins.species
  #    #  For each species the first and last bins must be non-empty, i.e.
  #    #  w_{s1}, w_{s,J_s +1} > 0. **Write code to check that before
  #    #  calling this function (since this gets repeatedly called).
  #    #dataBinForLikeSummary: table data frame with one row for each species,
  #    #  giving
  #    #  the minimum lower bound [w_{s1}] and maximum upper bound [w_{s,J_s +1}]
  #    #  and the number of counts for that species, where J_s is the number of
  #    #  bins for species s (won't need to explicilty specify). Columns are:
  #    #      SpecCode: code for each species [s]
  #    #      wminSpecies: minimum lower bound [w_{s1}]
  #    #      wmaxSpecies: maximum upper bound [w_{s,J_s +1}]
  #    #      n_s: total number of counts for species s [n_s]
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
  # **MOVE THESE TO PRE-PROCESS function:
  #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
  #       d[1] == 0 | d[J] == 0 | min(d) < 0)
  #       stop("Parameters out of bounds in negLL.PLB")
  if(b != -1)
      {  # From updated equation (A.63**[number will change]):

         temp2 = mutate(dataBinForLike,
                        comp2 = Number * log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
         comp2Sum = sum(temp2$comp2)

         logLL = - n * log( abs( xmax^(b+1) - xmin^(b+1) ) ) + comp2Sum
         neglogLL = - logLL      # Negative log-likelihood
      } else
      { stop("NOT DONE b=-1 yet; adapt, but first correct, negLL.PLB.binned")
      }
    return(neglogLL)
  }



