# PLBfunctions.r - functions to be sourced, including density, distribution
#  function and random number generator for unbounded and bounded power-law
#  distributions.
#
# CONTENTS
#
# Statistical functions:
# dPL - probability density function for unbounded power-law distribution
# pPL - probability distribution function for unbounded power-law distribution
# rPL - random numbers from an unbounded power-law distribution
# dPLB - probability density function for bounded power-law distribution
# pPLB - probability distribution function for bounded power-law distribution
# rPLB - random numbers from a bounded power-law distribution
# negLL.PLB - negative log-likelihood function for PLB model
# sum.bins - calculate the total sum of values within a bin
# Llin.method - fitting data using the Llin method
# LBNbiom.method - fitting data using the LBbiom and LBNbiom methods
# LBmizbinsFuns - calculate the bin breaks for the LBmiz method from mizer
# log2bins - construct bins that double in size and encompass the data
# eightMethods - computes exponents for a dataset using all eight methods, for
#  a second manuscript
# negLL.PLB.binned - negative log-likelihood function for PLB when data are
#  only available in binned form
# 
# Plotting functions:
# lm.line - plot straight line of lm fit but restricted to the x values
# gap.barplot.cust - customised version of Jim Lemon's gap.barplot for
#  histograms with a break in an axis
# qqtab - constructs automated LaTeX code for tables of quantiles
# confPlot - plotting of the confidence intervals for Figure 4
# histAxes - histogram axes for histogram plots of estimated b values
#  (Figure 3 and others)
# histAxes2 - histAxes adapted for fitting3rep-n10000.r, for n=10,000 sample size
# logTicks - add axes and tick marks to a log-log plot to represent
#  unlogged values (e.g. Figures 2(h) and 6(b))
# legJust - add legend to a plot
#
#  2nd Sept 2014 onwards.

require(plotrix)       # for axis.break function in gap.barplot.cust

# Statistical functions:

dPL = function(x = 1, b = -2, xmin = 1)    
  {
  # Computes probability density function for an unbounded power-law
  #  (Pareto) distribution
  #
  # Args:
  #   x: vector of values to compute the density function
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of probability density corresponding to vector x
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in dPL")
    C = - (b+1) / xmin^(b+1)
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = C * x[x >= xmin]^b
    return(y)                              
  }

pPL = function(x = 10, b = -2, xmin = 1)    
  {
  # Computes probability distribution function, P(X <= x),  for an
  #   unbounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in qPL")
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = 1 - (x[x >= xmin]/xmin)^(b+1)
    return(y)                              
  }

rPL = function(n = 1, b = -2, xmin = 1)    
  {
  # Computes random numbers from an unbounded power-law (Pareto) distribution
  #
  # Args:
  #   n: number of random numbers in sample. If 'length(n) > 1', the length is
  #       taken to be the number required. 
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of length n of independent random draws from the distribution
  #
  # Uses the inverse method (e.g. p1215 of Edwards 2008, Journal
  #   of Animal Ecology, 77:1212-1222).
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in rPL")
    u = runif(n)                            
    y = xmin * ( 1 - u ) ^ (1/(b+1))       
    return(y)                              
  }

# Bounded power-law distribution:

dPLB = function(x = 1, b = -2, xmin = 1, xmax = 100)    
  {
  # Computes probability density function for a bounded power-law
  #  (Pareto) distribution
  #
  # Args:
  #   x: vector of values to compute the density function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability density corresponding to vector x
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in dPLB")
    if(b != -1)
        { C = (b+1) / ( xmax^(b+1) - xmin^(b+1) )
        } else
        { C = 1/ ( log(xmax) - log(xmin) )
        }
    y = 0 * x     # so have zeros where x < xmin or x > xmax
    y[x >= xmin & x <= xmax] = C * x[x >= xmin & x <= xmax]^b
    return(y)                              
  }

pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100)    
  {
  # Computes probability distribution function, P(X <= x),  for a
  #   bounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x     # so have zeros where x < xmin
    y[x > xmax] = 1  # 1 for x > xmax
    if(b != -1)
        {  xmintobplus1 = xmin^(b+1)
           denom = xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] =
               ( x[x >= xmin & x <= xmax]^(b + 1) - xmintobplus1 ) / denom
        } else
        {  logxmin = log(xmin)
           denom = log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)                              
  }

rPLB = function(n = 1, b = -2, xmin = 1, xmax = 100)    
  {
  # Computes random numbers from a bounded power-law distribution
  #
  # Args:
  #   n: number of random numbers in sample. If 'length(n) > 1', the length is
  #       taken to be the number required. 
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  #
  # Returns:
  #   vector of length n of independent random draws from the distribution
  #
  # Uses the inverse method (e.g. p1215 of Edwards 2008, Journal
  #   of Animal Ecology, 77:1212-1222).
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in rPLB")
    u = runif(n)
    if(b != -1)
        { y = ( u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
        } else
        { y = xmax^u * xmin^(1-u)
        }
    return(y)                              
  }

negLL.PLB = function(b, x, n, xmin, xmax, sumlogx)
  {
  # Calculates the negative log-likelihood of the parameters b, xmin and xmax
  #  given data x for the PLB model. Returns the negative log-likelihood. Will
  #  be called by nlm or similar, but xmin and xmax are just estimated as the
  #  min and max of the data, not numerically using likelihood.
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   x: values of data x
  #   n: length(x), have as an input to avoid repeatedly calculating it 
  #   xmin: minimum value of x, have as an input to avoid repeatedly calculating
  #   xmax: maximum value of x, have as an input to avoid repeatedly calculating
  #   sumlogx: sum(log(x)) as an input, to avoid repeatedly calculating
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
  }

sum.bins = function(x, bin.breaks)
    {
    # Calculate the total sum of the values within each bin, e.g. total
    #  biomass of individuals within each size class. [May have used now used
    #  cut() instead for this].
    #
    # Args:
    #  x: vector of values (e.g. masses of individual fish)
    #  bin.breaks: breaks of the bins into which to partition the values,
    #   in ascending order
    # Returns:
    #  vector of totals of x within each bin
       N = length(bin.breaks) 
       if(min(x) < bin.breaks[1] | max(x) > bin.breaks[N])
           stop("Data out of range of bins in sum.bins")
       y = rep(NA, N-1)
       for(i in 1:(N-1))
           {  y[i] = sum( x[x >= bin.breaks[i] & x < bin.breaks[i+1]])
           }
       y[N-1] = y[N-1] + sum( x[ x == bin.breaks[i+1] ] )
           # since final bin needs to include max value
       return(y)
   }

Llin.method = function(bodyMass, num.bins = NULL, binBreaks = NULL)
    {
    # Use the Llin method, which is plotting binned counts on log-linear
    #  axes and then fitting a regression, as done by Daan et al. 2005.
    # 
    # Args:
    #  bodyMass: vector of individual body masses
    #  num.bins: number of bins to be used, though this is only a suggestion
    #   since can get over-ridden by hist().
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
        require(dplyr)   
        if(!is.vector(bodyMass)) stop("bodyMass not a vector in Llin.method")
        if(anyNA(bodyMass)) stop("bodyMass contains NA's in Llin.method")
        if(min(bodyMass) <= 0) stop("bodyMass needs to be >0 in Llin.method")
        x = bodyMass
        #
        if(!is.null(binBreaks))
           { 
            if(min(diff(binBreaks)) < 0) stop("binBreaks need to be increasing")
            breaks = binBreaks
           }  else { breaks = num.bins }
        hLlin = hist(x, breaks=breaks, plot=FALSE) # was breaks=num.bins, 4/11/15
        hLlin.mids = hLlin$mids

        hLlin.log.counts = log(hLlin$counts)
        hLlin.log.counts[ is.infinite(hLlin.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

        hLlin.lm = lm( hLlin.log.counts ~ hLlin.mids, na.action=na.omit)

        y = list(mids = hLlin.mids, log.counts = hLlin.log.counts,
            counts = hLlin$counts, lm = hLlin.lm, slope = hLlin.lm$coeff[2],
            breaks = hLlin$breaks,
            confVals = confint(hLlin.lm, "hLlin.mids",0.95))
        return(y)
       }

LBNbiom.method = function(bodyMass = NULL, counts = NULL,
       binBreaks = NULL, lowerCutOff = 0)
    {
    # Use the log-binning with normalisation technique (LBN method) to
    #  calculate the slope of the biomass size spectra. Slope is from fitting
    #  a linear regression of log10(normalised biomass in bin)
    #  against log10(midpoint of bin). Bins can be defined by user,
    #  else are created to double in size. Also calculates slope
    #  for biomasses not being normalised (LBbiom method).
    # 
    # Args:
    #  bodyMass: vector of individual body masses.
    #  counts: dataframe (or array) with first column being a body mass
    #   value, and second column being the counts of the number of individuals
    #   for that body mass. Only bodyMass or counts can be specified.
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
    #   
    #   indiv: dataframe with a row for each bodyMass value, where the
    #    columns are: x - original bodyMass values.
    #                 binMid, binMin, binMax, binWidth - midpoint, minimum,
    #                  maximum, and width, respectively, of the bin within
    #                  which the x value falls.
    #     If indiv has >=10^6 rows then it isn't saved.
    #     If counts was specified then, for now, an equivalent bodyMass
    #      vector was created and is column x (i.e. body masses are repeated).
    #
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
        if(!is.null(bodyMass) & !is.null(counts)) {
            stop("need only one of bodyMass or counts in LBNbiom.method") }
        if(is.null(bodyMass) & is.null(counts)) {
            stop("need bodyMass or counts in LBNbiom.method") }
        if(!is.null(bodyMass)) {
          if(!is.vector(bodyMass))stop("bodyMass not a vector in LBNbiom.method")
          if(anyNA(bodyMass)) stop("bodyMass contains NA's in LBNbiom.method")
          if(min(bodyMass) <= 0)stop("bodyMass needs to be >0 in LBNbiom.method")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in LBNbiom.method")
          if(min(counts[,1]) < 0) {
              stop("body masses in counts need to be >= 0 in LBNbiom.method") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in LBNbiom.method") }
          }
        # First wrote code for x = bodyMass, then wanted to add in the option
        #  to have counts as an input. Could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code.
        if(!is.null(bodyMass)) {x = bodyMass} else 
           {x = rep(counts[,1], counts[,2]) }
        #
        if(is.null(binBreaks))
           {
            binBreaks = 2^( floor(log2(min(x))) : ceiling(log2(max(x))) )
           } else
           { 
            if(min(diff(binBreaks)) < 0) { stop("binBreaks need to be increasing
                                              in LBNbiom.method")}
           }
        if(!(lowerCutOff %in% c(0, binBreaks))) { stop("need lowerCutOff
                           to be 0 or one of the binBreaks in LBNbiom.method") }
           # We could allow a user to specify a lowerCutoff that was not
           #  a binBreak, but it just makes defining the binBreaks a bit more
           #  fiddly -- code could be modified if a user wished. Although in
           #  practice people would plot the binned data and then choose which
           #  points (binned counts) to ignore when fitting the regression.
        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))        
        indiv$binMin =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-1])
        # indiv$binWidth =cut(x, breaks=binBreaks, right=FALSE,
        #    include.lowest=TRUE, labels = diff(binBreaks))
        # indiv = mutate(indiv, binWidth = binMax - binMin)
           # Above commands avoid any problems with bins with 0 counts.
           # Don't really need all of them, but include for completeness.
        indiv$binMid = as.numeric(as.character(indiv$binMid))
        indiv$binMin = as.numeric(as.character(indiv$binMin))
        indiv$binMax = as.numeric(as.character(indiv$binMax))
           # Now calculate biomass in each bin class:
        binVals = summarise(group_by(indiv, binMid),
            binMin = unique(binMin),
            binMax = unique(binMax), binWidth = binMax - binMin,
            totalBiom = sum(x), totalBiomNorm = totalBiom / binWidth )
           # binWidth uses new columns binMax and binMin
        binVals = binVals[order(binVals$binMid),]   # order by binMid
        #
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
        if(dim(indiv)[1] < 10^6) {       # only save indiv if not too big
          y = list(indiv = indiv, binVals = binVals,
            unNorm.lm = unNorm.lm, unNorm.slope = unNorm.slope,
            unNorm.conf = unNorm.conf,
            norm.lm = norm.lm, norm.slope = norm.slope, norm.conf = norm.conf,
            lowerCutOff = lowerCutOff)
          } else
          {    
          y = list(binVals = binVals,
            unNorm.lm = unNorm.lm, unNorm.slope = unNorm.slope,
            unNorm.conf = unNorm.conf,
            norm.lm = norm.lm, norm.slope = norm.slope, norm.conf = norm.conf,
            lowerCutOff = lowerCutOff)
          }
        return(y)
       }

LBmizbinsFun = function(beta, xmin, xmax, k)
    {
    # Calculates the bin breaks for the LBmiz method taken from mizer, given
    #  xmin and xmax (min and max of data) and the number of bins, k. To be
    #  minimised by nlm to calculate beta.
    #
    # beta: to be calculated,  log10 beta is the constant binwidth on log10
    #   scale. beta is solution to 0 = beta^(k-2) * (2 * beta-1) - xmax/xmin
    #                                = 2 * beta^(k-1) - beta^(k-2) - xmax/xmin
    # xmin: minimum of data (lower bound of lowest bin)
    # xmax: maximum of data (upper bound of highest bin)
    # k: number of bins
    # if( beta < 1) stop("beta needs to be >1; try reducing the number of
    #                       requested bins (k)")
    if( xmin <= 0 | xmin >= xmax | xmin >= xmax | k < 2 )
        { stop("Parameters out of bounds in LBmizbinsFun") }
    #fun = abs(beta^(k-2) * (2 * beta - 1) - xmax/xmin)   # to minimise
    # Use log scale for numerical stability
    fun = abs( (k-2)*log(beta) + log(2 * beta - 1) - log(xmax) + log(xmin))
                                         # to minimise
    return(fun)
    }


log2bins = function(x = NULL, counts = NULL)
    {
    #  Construct bins that double in size and encompass the data,
    #   resulting in binBreaks   ..., 0.25, 0.5, 1, 2, 4, 8, 16,....
    #   as necessary.
    #  Adapting from LBNbiom.method().
    # 
    # Args:
    #  x: vector of individual values (e.g. body masses).
    #  counts: dataframe (or array) with first column being an x value
    #  (e.g. body mass), and second column being the counts of the
    #   number of individuals for that value.
    #   Only x or counts can be specified.
    #
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
    #    columns are: binMid, binMin, binMax, binWidth - midpoint, minimum,
    #                  maximum, and width, respectively, of the bin
    #                 binCount - total number of individuals in that bin
    #                 binCountNorm - binCount / binWidth
    #                 binSum - sum of individual values in that bin
    #                  (appropriate if x represents biomass, but not length)
    #                 binSumNorm - binSum / binWidth    
    #                 log10.... - log10 of some of the above quantities
        require(dplyr)
        if(!is.null(x) & !is.null(counts)) {
            stop("need only one of x or counts in log2bins") }
        if(is.null(x) & is.null(counts)) {
            stop("need x or counts in log2bins") }
        if(!is.null(x)) {
          if(!is.vector(x))stop("x not a vector in log2bins")
          if(anyNA(x)) stop("x contains NA's in log2bins")
          if(min(x) <= 0)stop("x needs to be >0 in log2bins")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in log2bins")
          if(min(counts[,1]) < 0) {
              stop("x values in counts need to be >= 0 in log2bins") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in log2bins") }
          }
        # As for LBNbiom.method(), could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code.
        if(is.null(x))  
           {x = rep(counts[,1], counts[,2]) }
        #
        binBreaks = 2^( floor(log2(min(x))) : ceiling(log2(max(x))) )

        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))        
        indiv$binMin =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-1])
        # indiv$binWidth =cut(x, breaks=binBreaks, right=FALSE,
        #    include.lowest=TRUE, labels = diff(binBreaks))
        # indiv = mutate(indiv, binWidth = binMax - binMin)
           # Above commands avoid any problems with bins with 0 counts.
           # Don't really need all of them, but include for completeness.
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
        binVals = binVals[order(binVals$binMid),]   # order by binMid
        #
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

eightMethods = function(oneYear = 1980, 
  dataForYear = filter(data, Year == oneYear), figName = "eightMethods" )
  {
  # Computes exponents for a dataset using all eight methods. It was developed
  #  and is called in nSea15analysis.Snw, which was modified from what was
  #  in fitting2.r.  
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
  # Need x, a vector of individual fish sizes (lengths or weights), which is how
  #  the original methods functions are written.
  # Now adding explicit methods in eightMethods.counts, such as
  #  Llin.method.counts, to deal explicitly with counts, and that
  #  should also work for non-integer counts. For integer
  #  counts, the expansion to give x should give the same results. 
  x = rep(dataForYear$bodyMass, dataForYear$Number)
              # 3.3 million for 1980 for old nSea15, but it was quick
  # return(head(dataForYear))   }       # when testing
  log.x = log(x)                     
  sum.log.x = sum( log.x ) 
  xmin = min(x)
  xmax = max(x)

  figheight = 7 # 5.6     For 4x2 figure
  figwidth = 5.7    # 5.7 inches for JAE

  num.bins = 8   # number of bins for standard histogram and Llin method, though
                 #  this is only a suggestion (and can get overridden). Daan used
                 #  8 bins.

  # postscript("nSea1980-fitting2a.eps", height = figheight,
  #           width = figwidth, horizontal=FALSE, paper="special")
  # Postscript plots all 3million points, ends up >150Mb file. So try .png:
  png(paste(figName, "-", oneYear, ".png", sep=""), height = figheight,
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
  hLlin.list = Llin.method(x, num.bins = num.bins)

  #eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
  #    "b" = hLlin.list$slope,
  #    "confMin"= hLlin.list$confVals[1], "confMax"= hLlin.list$confVals[2]) 
  # That gives hLlin.mids as the name of the row, so does this though
  hLlin.b = hLlin.list$slope
  hLlin.confMin = hLlin.list$confVals[1]
  hLlin.confMax = hLlin.list$confVals[2]
  eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
      "b" = hLlin.b, "confMin"= hLlin.confMin, "confMax"= hLlin.confMax,
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

  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LT"),
      "b" = hLT.slope, "confMin"= hLT.conf[1], "confMax"= hLT.conf[2], 
      row.names=NULL))
   
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

  hLTplus1.conf = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95)
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LTplus1"),
      "b" = hLTplus1.slope, "confMin"= hLTplus1.conf[1], 
      "confMax"= hLTplus1.conf[2], row.names=NULL))
  
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
  
  hLBmiz.num.bins = num.bins
  
  beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est
  
  # hLBmiz.bins = c(beta^(0:(k-1)) * xmin, xmax)
  hLBmiz.bins = c(beta^(0:(hLBmiz.num.bins-1)) * min(x), max(x))
   # Mizer bin specification, with final bin being same width as penultimate bin
  hLBmiz = hist(x, breaks=hLBmiz.bins, plot=FALSE)     # linear scale. Only for counts.

  # From mizer's getCommunitySlopeCode.r:
  #  "Calculates the slope of the community abundance through time by
  #  performing a linear regression on the logged total numerical abundance
  #  at weight and logged weights (natural logs, not log to base 10, 
  #  are used)."  So regress log(total counts) against log(weights)
  #  (not log10 and not normalised). And it's actually
  #   on the minima of the bins (their w).

  hLBmiz.log.min.of.bins = log(hLBmiz.bins[-length(hLBmiz.bins)])# min of bins
  hLBmiz.log.counts = log(hLBmiz$counts)
  hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA  
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)

  # hLBmiz.slope = hLBmiz.lm$coeff[2]
  # Need to subtract 1, since want to work in terms of b not slopes now
  hLBmiz.b = hLBmiz.lm$coeff[2] - 1

  hLBmiz.conf = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95) - 1
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBmiz"),
      "b" = hLBmiz.b, "confMin"= hLBmiz.conf[1], 
      "confMax"= hLBmiz.conf[2], row.names=NULL))

  plot( hLBmiz.log.min.of.bins, hLBmiz.log.counts,
     xlab=expression(paste("Log (minima of bins for data ", italic(x), ")")),
     ylab = "Log (count)", mgp=mgpVals)
     # axes=FALSE, xaxs="i", yaxs="i", , xlim=c(log10(1), log10(650)), 
     #  ylim=c(log10(0.7), log10(1100)))  # So axes are logged

  lm.line(hLBmiz.log.min.of.bins, hLBmiz.lm)
  legend("bottomleft", paste("(d) LBmiz b=", signif(hLBmiz.b, 3)), 
         bty="n", inset=inset)

  # mizer biomass size spectra - see option (not using) in mizerBiom.eps below.


  # LBbiom method - binning data using log2 bins, calculating biomass not counts
  #  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
  #  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

  hLBNbiom.list = LBNbiom.method(x)    # Does this method and the next.

  hLBbiom.b = hLBNbiom.list[["unNorm.slope"]] - 2
  hLBbiom.conf = hLBNbiom.list[["unNorm.conf"]] - 2
  
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBbiom"),
      "b" = hLBbiom.b, "confMin"= hLBbiom.conf[1], 
      "confMax"= hLBbiom.conf[2], row.names=NULL))
  
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

  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LBNbiom"),
      "b" = hLBNbiom.b, "confMin"= hLBNbiom.conf[1], 
      "confMax"= hLBNbiom.conf[2], row.names=NULL))
  
  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])
  legend("bottomleft", paste("(f) LBNbiom b=",
     signif(hLBNbiom.b, 3)), bty="n", inset=inset)

  # Cumulative Distribution, LCD method
  x.sorted = sort(x, decreasing=TRUE)
  logSorted = log(x.sorted)
  logProp = log((1:length(x))/length(x))

  hLCD.lm = lm(logProp ~ logSorted)   # plot(fitsortedlog10) shows
                                                   #  residuals not good
  hLCD.slope = hLCD.lm$coeff[2]

  hLCD.b = hLCD.lm$coeff[2] - 1
  hLCD.conf = confint(hLCD.lm, "logSorted", 0.95) - 1

  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("LCD"),
      "b" = hLCD.b, "confMin"= hLCD.conf[1], 
      "confMax"= hLCD.conf[2], row.names=NULL))
    
  plot(logSorted, logProp, main="",
     xlab=expression(paste("Log ", italic(x))),
     ylab=expression( paste("Log (prop. of ", values > italic(x), ")")),
     mgp=mgpVals) # , axes=FALSE)
     #xlim=c(0.8, 1000), xaxs="i", ylim=c(0.0008, 1), yaxs="i",

  lm.line(logSorted, hLCD.lm, col="red")
  # murankfreq = 1 - fitsortedlog10$coeff[2]       # mu = 1 - slope
  legend("bottomleft", paste("(g) LCD b=", signif(hLCD.b, 3)), bty="n",
       inset=inset)
 
  # MLE (maximum likelihood method) calculations.
  
  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1
      
  PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
      xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )
  
  PLB.bMLE = PLB.minLL$estimate
  
  # 95% confidence intervals for MLE method.
  
  PLB.minNegLL = PLB.minLL$minimum
  
  # Values of b to test to obtain confidence interval. For the real movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.
  
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.001)  # If make 0.0001 then do 
                              # get an interval for raw 1980 data 
      
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
          PLB.LLvals[i] = negLL.PLB(bvec[i], x=x, n=length(x), xmin=xmin,
              xmax=xmax, sumlogx=sum.log.x)   
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
  
  # MLE.rep.xmax[iii] = xmax   - not storing xmax for now
  eightMethodsRes = rbind(eightMethodsRes, 
      data.frame("Year"=oneYear, "Method"=as.factor("MLE"),
      "b" = PLB.bMLE, "confMin"= PLB.MLE.bConf[1], 
      "confMax"= PLB.MLE.bConf[2], row.names=NULL))
  
  # To plot rank/frequency style plot:
  plot(sort(x, decreasing=TRUE), 1:length(x), log="xy",
       xlab=expression(paste("Values, ", italic(x))),
       ylab=expression( paste("Number of ", values >= x)), mgp=mgpVals)
    # , xlim=c(2, 500), ylim=c(0.8, 40000), axes=FALSE, xaxs="i", yaxs="i", 
    #  mgp=mgpVals)

  x.PLB = seq(min(x), max(x), length=1000)     # x values to plot PLB
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                    xmax = max(x.PLB))) * length(x)    
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

negLL.PLB.binned = function(b, w, d, J=length(d), xmin=min(w), xmax=max(w))
                                        # sumlogx)
  {
  # Calculates the negative log-likelihood of the parameters b, xmin and xmax
  #  given binned data for the PLB model. Returns the negative log-likelihood.
  #  Will be called by nlm or similar, but xmin and xmax will just be estimated
  #  as the min of lowest bin and max of the largest bin (that is their MLEs),
  #  no need to do numerically. See Appendix for derivation.
  #      
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   w: vector of length J+1 giving the bin breaks w_1, w_2, ..., w_{J+1}
  #   d: vector of length J giving the count in each bin. Must have d_1, d_J >0
  #   J: number of bins (length of d)
  #   xmin: minimum value of bins, as an input to avoid repeatedly calculating
  #   xmax: maximum value of bins, as an input to avoid repeatedly calculating
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
         d[1] == 0 | d[J] == 0 | min(d) < 0)
         stop("Parameters out of bounds in negLL.PLB")
    n = sum(d)
    if(b != -1)
      { neglogLL = n * log( abs( w[J+1]^(b+1) - w[1]^(b+1) ) ) -
            sum( d * log( abs( w[-1]^(b+1) - w[-(J+1)]^(b+1) ) ) ) 
      } else
      { neglogLL = 1 / (log(w[J+1]) - log(w[1]) ) *
            sum( d * ( log(w[-1]) - log(w[-(J+1)]) ) )
      }
    return(neglogLL)
  }


# Plotting functions:

lm.line = function(x.vector, lm.results, ...)
  {
  # Plots a straight line between lowest and highest values of x.vector
  #  based on lm results, to use instead of abline(lm.output)
  #  which extends beyond the fitted data.
  #
  # Args:
  #   x.vector: values of x (that have been fitted to)
  #   lm.results: output from an lm fit, with intercept lm.results$coeff[1]
  #     and gradientslope lm.results$coeff[2]
  #   ...: arguments to lines
  # Returns:
  #   adds a line to an existing plot
     intercept = lm.results$coeff[1]
     grad = lm.results$coeff[2]
     lines(c( min(x.vector), max(x.vector)), c(grad * min(x.vector) + intercept, 
      grad * max(x.vector) + intercept), ...)
  }

gap.barplot.cust = 
function (y, gap, xaxlab, xtics, yaxlab, ytics, midpoints, breakpoints,
    ylim = NA, xlab = NULL, 
    ylab = NULL, horiz = FALSE, col = NULL, xlim, N = n, ...) 
    {
    # For Figure 1, for a histogram (barplot) with a gap in the y-axis.
    #  Customising gap.barplot from the package plotrix by Jim Lemon.
    #  Several options here are customised for the particular plot (and to
    #  change a few of the defaults in gap.barplot) so the code
    #  requires some modifiying to use more generally. 
    if (missing(y)) 
        stop("y values required")
    # x <- 1:length(y)        # Original
    x <- midpoints            # AE adding; midpoints is defined in raw1infexamp.r
    if (missing(gap)) 
        stop("gap must be specified")
    if (is.null(ylab)) 
        ylab <- deparse(substitute(y))
    if (is.null(col)) 
        col <- color.gradient(c(0, 1), c(0, 1, 0), c(1, 0), length(y))
    else if (length(col) < length(y)) 
        rep(col, length.out = length(y))
    littleones <- which(y <= gap[1])
    bigones <- which(y >= gap[2])
    if (any(y > gap[1] & y < gap[2])) 
        warning("gap includes some values of y")
    gapsize <- gap[2] - gap[1]
    if (missing(xaxlab)) 
        xaxlab <- as.character(x)
    # xlim <- c(min(x) - 0.4, max(x) + 0.4)    # Original
    xlim <- xlim
    if (is.na(ylim[1])) 
        ylim <- c(min(y), max(y) - gapsize)
    if (missing(ytics)) 
        ytics <- pretty(y)
    if (missing(yaxlab)) 
        yaxlab <- ytics
    littletics <- which(ytics < gap[1])
    bigtics <- which(ytics >= gap[2])
    halfwidth <- min(diff(x))/2
    if (horiz) {                   # AE not editing anything for horizontal
        if (!is.null(xlab)) {
            tmplab <- xlab
            xlab <- ylab
            ylab <- tmplab
        }
        plot(0, xlim = ylim, ylim = xlim, ylab = ylab, axes = FALSE, 
            type = "n", ...)
        plot.lim <- par("usr")
        botgap <- ifelse(gap[1] < 0, gap[1], plot.lim[1])
        box()
        axis(2, at = x, labels = xaxlab)
        axis(1, at = c(ytics[littletics], ytics[bigtics] - gapsize), 
            labels = c(ytics[littletics], ytics[bigtics]))
        rect(botgap, x[y < gap[1]] - halfwidth, y[y < gap[1]], 
            x[y < gap[1]] + halfwidth, col = col[y < gap[1]])
        rect(botgap, x[bigones] - halfwidth, y[bigones] - gapsize, 
            x[bigones] + halfwidth, col = col[bigones])
        axis.break(1, gap[1], style = "gap")
    }
    else {                        # AE editing here for vertical bars
        plot(0, xlim = xlim, ylim = ylim, ylab = ylab, axes = FALSE, xlab=xlab,
            type = "n", ...)               # AE added xlab=xlab
        plot.lim <- par("usr")
        botgap <- ifelse(gap[1] < 0, gap[1], plot.lim[3])
        box()
        # axis(1, at = x, labels = xaxlab)  # - original was x=1,2,3,...,58, is
                                            #  now the midpoints, but I want
                                            #  ticks at breakpoints,
                                            #  hmadeup$breaks
        axis(1, breakpoints, tcl=-0.2, labels=rep("", length(breakpoints)))
                                                # short ones at breaks
        # axis(1, at=seq(0, 600, 50), labels = rep("", 13), mgp=c(1.7,0.6,0))
                                             # long ticks where I want
        axis(1, at=seq(0, 1000, 100),
             labels = seq(0, 1000, 100), mgp=c(1.7,0.6,0))
                  # long ticks labelled where I want, not automated though
        axis(2, at = c(ytics[littletics], ytics[bigtics] - gapsize), 
            labels = c(ytics[littletics], ytics[bigtics]), mgp=c(1.7,0.6,0))
                       # AE adding mgp, these are the ones with numbers on
        # Short tics:
        bottomShortTics = seq(0, gap[1], 2)        # below gap
        axis(2, bottomShortTics, tcl=-0.2,
             labels=rep("", length(bottomShortTics)))    
        topShortTics = seq(gap[2], N, 2)-(gap[2]-gap[1])  # above gap
        axis(2, topShortTics, tcl=-0.2,
             labels=rep("", length(topShortTics)))
        rect(x[y < gap[1]] - halfwidth, botgap, x[y < gap[1]] + 
            halfwidth, y[y < gap[1]], col = col[y < gap[1]])
        rect(x[bigones] - halfwidth, botgap, x[bigones] + halfwidth, 
            y[bigones] - gapsize, col = col[bigones])
        axis.break(2, gap[1], style = "gap")
    }
    invisible(x)
}

qqtab = function(xx, dig=2, true=b.known, quants = c(0.25, 0.75))   
  {
    # Quantile table function, to construct lines of quantiles to copy
    #  into LaTeX. Does quants[1], 50% (median), mean, quants[2] quantiles,
    #  and the %age of values < the true value.
    # xx: vector of values to give quantiles for.
    # dig: number of decimal places to give.
    # true: the true value of the quantity being estimated, will
    #  depend on method.
    # quants: quantiles to use, with defaults of 25% and 75%.
      paste( c( prettyNum(round(quantile(xx, quants[1]), 
                                  digits=dig), big.mark=","), 
         " & ", prettyNum(round(quantile(xx, 0.50), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(mean(xx), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(quantile(xx, quants[2]), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(sum(xx < true)/100, digits=0),
                                  big.mark=",")), 
         sep="", collapse="")             # , "\\%"),
  }  




confPlot = function(repConf, legName, b.true = b.known, inCol="darkgrey",
    outCol="blue", vertCol="red", pchVal = 20, cexVal = 0.0, xLim = NULL,
    colourCode = TRUE, vertThick = 1, horizThick = 0.3,
    thin = 33, horizLines = FALSE,
    horizLinesOut = TRUE, horizLinesIn = TRUE, yLab = "Sample number",
    yTicks = seq(0, 300, 50), yLabels = TRUE, vertFirst = FALSE,
    insetVal = c(-0.08, -0.06), insetVal2 = c(-0.08, 0.07), xsmallticks=NULL)
    {
    # Plotting function for confidence intervals of the repeated estimates
    #  for one method. Gets called eight times to produce Figure 4. Plots
    #  horizontal lines for the intervals, colour coded as to whether they
    #  include the true value of b.
    # Args:
    #  repConf: two-column data-frame, columns are confMin and confMax, the
    #   minimum and maximum of the 95% confidence interval for b,
    #   and rows correspond to each simulated data set. When
    #   confPlot is called, for some methods b = slope - 1 or slope -2, which
    #   gets done in the call. Originally plotted the slopes, but converting
    #   to b makes more sense.
    #  legName: legend name for that panel
    #  b.true: what b should be
    #  inCol: colour for intervals that b.true is inside
    #  outCol: colour for intervals that b.true is outside
    #  vertCol: colour of vertical line for b.true
    #  pchVal: pch for points at endpoints of intervals
    #  cexVal: size of points (default of 0 does not plot them)
    #  xLim: xlim range
    #  colourCode: colour code the figure (may not completely work)
    #  vertThick: thickness of vertical line for b.true
    #  horizThick: thickness of horizontal lines vertical line for b.true
    #  thin: number of values to thin the values for plotting. Only works for
    #   33 or 99 for now (since they are factors of 9999).
    #  horizLines: whether to plot horizontal grey lines or not as example
    #     (not needed now since doing lines for intervals)
    #  horizLinesOut, horizLinesIn: whether to plot horizontal lines for
    #     intervals for which true value is outside/inside the interval.
    #  yLab: label for y axis
    #  yTicks: where to have tickmarks on y-axis
    #  yLabels: whether or not to label tickmarks on y-axis
    #  vertfirst: whether or not to plot vertical line first in order to see
    #   the horizontal lines better (for LCD plot)
    #  insetVal: inset shift for naming the panel
    #  insetVal2: inset shift for printing observed coverage percentage
    #  xsmallticks: where to put unlabelled small tick marks on x-axis
    if(!colourCode) outCol = inCol
    if(!(thin %in% c(33,99))) stop("Need to edit confPlot if thin not 33 or 99")
    if(is.null(xLim))
        {
            rangeVal = range(repConf)
            xLim = c(floor(rangeVal[1]), ceiling(rangeVal[2]))
                  # min(repConf[,1])), ceiling(max(repConf[,2])))
        }
    repConf = mutate(repConf,
      inConf = (b.true > confMin & b.true < confMax))  # does CI cover b.true?

    repConf = mutate(repConf, confCol = NA)
    # Couldn't figure this out in dplyr:
    repConf[which(repConf$inConf), "confCol"] = inCol
    repConf[which(!repConf$inConf), "confCol"] = outCol
    repConf[, "num.rep"] = as.numeric(row.names(repConf))
              # To preserve the iteration number in case needed later

    repConf.sort.full = arrange(repConf, confMin)   # arranged by min of CI
    sum.inConf = sum(repConf.sort.full$inConf, na.rm=TRUE) /
        dim(repConf.sort.full)[1]      # get NA's for Llin, LT with xmax=10000
    # subset for better plotting:
    if(thin == 33)
        { repConf.sort = repConf.sort.full[c(1, 1+thin*(1:303)), ] }else
        { repConf.sort = repConf.sort.full[c(1, 1+thin*(1:101)), ] }
    # print(head(repConf.sort))
    # repConf.sort[, "num.sorted"] = as.numeric(row.names(repConf.sort))
    #  That doesn't work since it preserves the repConf.sort.full row names
    repConf.sort[, "num.sorted"] = 1:(dim(repConf.sort)[1])
              # To preserve the new row number since needed later
    plot(repConf.sort$confMin, repConf.sort$num.sorted, col=repConf.sort$confCol,
     xlim = xLim, ylim = c(0, 1.02*dim(repConf.sort)[1]),
     pch=pchVal, cex=cexVal, 
     xlab = "", 
     ylab = yLab, yaxt="n")
    if(vertFirst) {abline(v=b.true, lwd=vertThick, col=vertCol)}
    axis(2, at = yTicks, labels = yLabels, tck=-0.04)
    if(!is.null(xsmallticks))
        { axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)}
    # xlab = expression(paste("Estimate of slope (or ", italic(b), ")")),
    points(repConf.sort$confMax, repConf.sort$num.sorted,
           col=repConf.sort$confCol, pch=pchVal, cex=cexVal)
    legend("topleft", legName, bty="n", inset=insetVal)
    legend("topleft", paste(round(sum.inConf*100, dig=0), "%", sep=""),
           bty="n", inset=insetVal2)
    # legend("topleft", "hello", bty="n", inset=c(-0.08, -0.2))
    
    # Plot just the rows for which true value lie outside CI
    if(horizLinesOut)
        {
        outConf = filter(repConf.sort, (!inConf))    
        for(jj in 1:(dim(outConf)[1]))
          {
        
            lines(c(outConf$confMin[jj], outConf$confMax[jj]),
              c(outConf$num.sorted[jj], outConf$num.sorted[jj]),
              col=outCol, lwd=horizThick)
          }
        }

    # Plot just the rows for which true value lie inside CI
    if(horizLinesIn)
        {
        in.Conf = filter(repConf.sort, (inConf))    
        for(jj in 1:(dim(in.Conf)[1]))
          {
        
            lines(c(in.Conf$confMin[jj], in.Conf$confMax[jj]),
              c(in.Conf$num.sorted[jj], in.Conf$num.sorted[jj]),
              col=inCol, lwd=horizThick)
          }
        }
    
    # Plot equally spaced horizontal lines between 2.5 and 97.5 values
    #  (though even from the likelihood ratio test for MLE a 95% CI
    #  isn't actually the 2.5-97.5 range, it's a 95% interval).
    if(horizLines)
        {
        # row numbers to use to plot equally spaced horizontal lines
        horLineInd = round(c(2.5, 21.5, 40.5, 59.5, 78.5, 97.5) * 0.01 *
            dim(repConf.sort)[1])
        horiz.lines = repConf.sort[horLineInd,] 
        for(jj in 1:length(horLineInd))
            {
              lines(c(horiz.lines$confMin[jj], horiz.lines$confMax[jj]),
              c(horiz.lines$num.sorted[jj], horiz.lines$num.sorted[jj] ),
              col="grey", lwd=horizThick)
            }
        }
    if(!vertFirst) {abline(v=b.true, lwd=vertThick, col=vertCol)}
    return(repConf.sort.full)     # repConf.sort is what's plotted though
}

histAxes = function()
    # Do the histogram axes in fitting3rep.r and later, since almost all
    #  panels will have same axes. Not very flexible, use histAxes2()
    #  to plot up to 10,000.
    {
    axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks
    axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)
    axis(2, at=c(0, 2000, 4000),
         labels = c(0, 2000, 4000),     
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
      axis(2, at=seq(0, 6500, 1000),
         labels = rep("", 7),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
      axis(2, at=seq(0, 6500, 500),
         labels = rep("", 14), mgp=c(1.7,0.7,0), tcl=-0.2)  # small ticks
      abline(v=b.known, col=vertCol, lwd=vertThick) 
}

histAxes2 = function()
    # Do the histogram axes for, in particular, fitting3rep-n10000.r,
    #  since the larger n sample size means affects the resulting histograms
    #  of estimated b. Not very flexible, just doing for the one figure.
    #  Copying what is used for Llin method (so could go back and use this
    #  function throughout earlier code).
    {
      axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0),
           cex.axis=cexAxis)  # big ticks
      axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)),
           tcl=-0.2)
      axis(2, at=c(0, 5000, 10000),
         labels = c(0, 5000, 10000),
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
      axis(2, at=seq(0, 10000, 1000),
         labels = rep("", 11),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
      abline(v=b.known, col=vertCol, lwd=vertThick)
    }


logTicks = function(xLim, yLim = NULL, tclSmall = -0.2, xLabelSmall = NULL,
       yLabelSmall = NULL, xLabelBig = NULL, mgpVal=c(1.6,0.5,0))
    {
    # Add axes and tick marks to a log-log plot to represent unlogged values.
    # Args:
    #  xLim: the x limits for the plot (unlogged scale); if NULL then
    #         do not add anything to x-axis
    #  yLim: the y limits for the plot (unlogged scale); if NULL then
    #         do not add anything to y-axis
    #  tclSmall: size of small tick marks
    #  xLabelSmall: which small tick marks on x-axis to label
    #  yLabelSmall: which small tick marks on y-axis to label
    #  xLabelBig: which big tick marks on the x-axis to label
    #   (when automated they can overlap, so may need to specify).
    #  mgpVal: mgp values for axes. See ?par 
    # Returns:
    #  Adds axes and big and small tick marks to the plot. Returns NULL.
    # Example:
    #  Adapt the following:
    #  plot(..., log="xy", xlab=..., ylab=..., xlim=..., ylim=..., axes=FALSE)
    #  xLim = 10^par("usr")[1:2]
    #  yLim = 10^par("usr")[3:4]
    #  logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500)) 
    # 
    ll = 1:9
    log10ll = log10(ll)
    box()
    # x axis
    if(!is.null(xLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      xEncompassLog = c(floor(log10(xLim[1])), ceiling(log10(xLim[2])))
      xBig = 10^c(xEncompassLog[1]:xEncompassLog[2])
      # Big unlabelled, always want these:
      axis(1, at= xBig, labels = rep("", length(xBig)), mgp = mgpVal)
      # Big labelled:
      if(is.null(xLabelBig)) { xLabelBig = xBig }
      axis(1, at= xLabelBig, labels = xLabelBig, mgp = mgpVal)
      # axis(1, at=c(1, 10, 100), labels = c(1, 10, 100), mgp=c(1.7,0.7,0))
      # Small unlabelled:
      axis(1, xBig %x% ll, labels=rep("", length(xBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(xLabelSmall))
          {
          axis(1, at=xLabelSmall, labels=xLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }
    # Repeat for y axis:
    if(!is.null(yLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      yEncompassLog = c(floor(log10(yLim[1])), ceiling(log10(yLim[2])))
      yBig = 10^c(yEncompassLog[1]:yEncompassLog[2])
      # Big labelled:
      axis(2, at= yBig, labels = yBig, mgp = mgpVal)
      # Small unlabelled:
      axis(2, yBig %x% ll, labels=rep("", length(yBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(yLabelSmall))
          {
          axis(2, at=yLabelSmall, labels=yLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }     
}

legJust = function(textVec, pos="topright", textWidth = "slope=-*.**",
        inset=0, logxy=FALSE)
    {
    # Add legend with right-justification, functionalising Uwe Ligges'
    #  example in ?legend. Really a way of adding text automatically in
    #  the corner (which legend() works out the positioning for).
    # Args:
    #  textVec: text for the legend, one element for each row
    #  pos:  position of the legend (not tested for all positions)
    #  inset: inset distance
    #  textWidth: width to make the text
    #  logxy: TRUE if axes are logarithmic
    # Returns:
    #  Adds legend to exiting plot. Returns NULL.
    # Example:
    #  plot(10:1)
    #  legJust(c("Method", paste("value=", mean(1:10), sep=""), "x=7.0"))
    leg = legend(pos, legend = rep(" ", length(textVec)),
               text.width = strwidth(textWidth), bty="n", inset=inset)
    if(logxy)
        { text(10^(leg$rect$left + leg$rect$w), 10^(leg$text$y), textVec, pos = 2) } else
        { text(leg$rect$left + leg$rect$w, leg$text$y, textVec, pos = 2) }
}

