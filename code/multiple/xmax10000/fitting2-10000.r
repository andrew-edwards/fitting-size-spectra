# fitting2-10000.r - redoing simulations with xmax=10000. 17th July 2015.
# fitting2.r - replacing the fitting methods with functions, that are all
#  in PLBfunctions.r. 5th June 2015.
# fitting1.r - simulating data then fitting using various techniques.
#  Somewhat adapated from raw1infexamp.r from JAE paper. 3rd September 2014.
#  Also used some of codeAbigeye1.r from PLOS ONE paper. 

# raw1infexamp.r - just one simulated data set, to do as example of each
#  fitting technique for my Figure 1. 26th June 2007
# raw1infpap.r - doing again for paper, then saving properly and making figures#  in raw1infpapfig.r. Also adding the cumulative distribution technique, as
#  Sims did also try that. 11th May 2007 Taken some from huntersinf2.r, mainly
#  CI for pow law inf, so some stuff has inf appended but some doesn't. In
#  huntersinf2.r if there was no inf it was assumed bounded. Maybe worth
#  rectifying that here.
# raw1inf.r - doing for unbounded power law. Like I had a while ago then.
#  24th April 2007 ....

rm(list=ls())
print(date())
source("../../PLBfunctions.r")  # to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy)

n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 10000         # known fixed value of xmax

# Sample from PLB distribution:
set.seed(42)      # To get the same observations for each run of code.
                  # 8 bins, only up to 400 with 2 empty.
newdata = 1      # 1 if generate new data, else 0 if load in previous set
                 # Do NOT change without doing new filenames
if(newdata == 1)
  {
  x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)
  } else
  {
  load(file="raw1infexamp10000.RData")    # or load in data for x
  }

# x is a vector of individual fish sizes (lengths or weights)

log.x = log(x)                      # to avoid keep calculating
sum.log.x = sum( log.x ) 
xmin = min(x)
xmax = max(x)

figheight = 7 # 5.6     For 4x2 figure
figwidth = 5.7    # 5.7 inches for JAE

num.bins = 8    # number of bins for standard histogram and Llin method, though
                #  this is only a suggestion (and can get overridden). Daan used
                #  8 bins.

postscript("fitting2a10000.eps", height = figheight, width = figwidth, horizontal=FALSE,         paper="special")  #  height=4, width=4   was 7,4
par(mfrow=c(4,2))
oldmai = par("mai")    #  0.95625 0.76875 0.76875 0.39375  inches I think,
                       #   think may be indpt of fig size
par(mai=c(0.4, 0.5, 0.05, 0.3))  # Affects all figures if don't change again
mgpVals = c(1.6,0.5,0)            # mgp values   2.0, 0.5, 0

# Notation:
# hAAA - h(istrogram) for method AAA.

# Llin method - plotting binned data on log-linear axes then fitting regression,
#  as done by Daan et al. 2005.
hLlin.list = Llin.method(x, num.bins = num.bins)

plot( hLlin.list$mids, hLlin.list$log.counts, xlim=c(0, 100*ceiling(xmax/100)),
     xlab=expression(paste("Bin midpoints for data ", italic(x))),
     ylab = "Log (count)", mgp=mgpVals)

lm.line(hLlin.list$mids, hLlin.list$lm)
inset = c(0, -0.04)     # inset distance of legend
legend("topright", paste("(a) Llin slope=", signif(hLlin.list$slope, 3)),
       bty="n", inset=inset)

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

plot( hLT.log.mids, hLT.log.counts, xlab=expression(paste("Log (bin midpoints for data ", italic(x), ")")), ylab = "Log (count)", mgp=mgpVals)

lm.line(hLT.log.mids, hLT.lm)
legend("topright", paste("(b) LT slope=", signif(hLT.slope, 3)), bty="n",
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
hLTplus1.lm = lm( hLTplus1.log10.counts ~ hLTplus1.log10.mids, na.action=na.omit)
hLTplus1.slope = hLTplus1.lm$coeff[2]

plot( hLTplus1.log10.mids, hLTplus1.log10.counts, xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")), ylab = "Log10 (count+1)", mgp=mgpVals)
 
lm.line(hLTplus1.log10.mids, hLTplus1.lm)
legend("topright", paste("(c) LTplus1 slope=", signif(hLTplus1.slope, 3)),
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

hLBmiz = hist(x, breaks=hLBmiz.bins, plot=FALSE)     # linear scale

# From mizer's getCommunitySlopeCode.r:
#  "Calculates the slope of the community abundance through time by performing a linear regression on the logged total numerical abundance at weight and logged weights (natural logs, not log to base 10, are used)."  So regress log(total
#  counts) against log(weights) (not log10 and not normalised). And it's actually
#  on the minima of the bins (their w).

hLBmiz.log.min.of.bins = log(hLBmiz.bins[-length(hLBmiz.bins)]) # min of each bin
hLBmiz.log.counts = log(hLBmiz$counts)
hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA  
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)
# Checking result for using midpoints rather than mins of bins, get slope -1.130:
# hLBmiz.lm = lm( hLBmiz.log.counts ~ log(hLBmiz$mids), na.action=na.omit)
hLBmiz.slope = hLBmiz.lm$coeff[2]

plot( hLBmiz.log.min.of.bins, hLBmiz.log.counts,
     xlab=expression(paste("Log (minima of bins for data ", italic(x), ")")),
     ylab = "Log (count)", mgp=mgpVals)
     # axes=FALSE, xaxs="i", yaxs="i", , xlim=c(log10(1), log10(650)), ylim=c(log10(0.7), log10(1100)))  # So axes are logged

lm.line(hLBmiz.log.min.of.bins, hLBmiz.lm)
legend("topright", paste("(d) LBmiz slope=", signif(hLBmiz.slope, 3)), bty="n",
       inset=inset)

# mizer biomass size spectra - see option (not using) in mizerBiom.eps below.


# LBbiom method - binning data using log2 bins, calculating biomass (not counts)
#  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
#  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

hLBNbiom.list = LBNbiom.method(x)    # Does this method and the next.

plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiom,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (biomass)", mgp=mgpVals)

lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["unNorm.lm"]])
legend("bottomleft", paste("(e) LBbiom slope=",
     signif(hLBNbiom.list[["unNorm.slope"]], 3)), bty="n", inset=c(-0.08, -0.04))

# LBNbiom method - on biomass, not counts, as per Julia's 2005 paper.
#  log2 bins of bodymass, sum the total biomass in each bin, normalise
#  biomasses by binwidths, fit regression to log10(normalised biomass) v
#  log10(midpoint of bin).

# hLBNbiom.list = LBNbiom.method(x) - already done above

plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals)

lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])
legend("topright", paste("(f) LBNbiom slope=",
     signif(hLBNbiom.list[["norm.slope"]], 3)), bty="n", inset=inset)

# Cumulative Distribution, LCD method
x.sorted = sort(x, decreasing=TRUE)
logSorted = log(x.sorted)
logProp = log((1:length(x))/length(x))

plot(logSorted, logProp, main="",
     xlab=expression(paste("Log ", italic(x))),
     ylab=expression( paste("Log (prop. of ", values >= italic(x), ")")),
     mgp=mgpVals) # , axes=FALSE)
     #xlim=c(0.8, 1000), xaxs="i", ylim=c(0.0008, 1), yaxs="i",

hLCD.lm = lm(logProp ~ logSorted)   # plot(fitsortedlog10) shows
                                                 #  residuals not good
hLCD.slope = hLCD.lm$coeff[2]
lm.line(logSorted, hLCD.lm, col="red")
# murankfreq = 1 - fitsortedlog10$coeff[2]       # mu = 1 - slope
legend("topright", paste("(g) LCD slope=", signif(hLCD.slope, 3)), bty="n",
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

bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.001) 
    
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

# To plot rank/frequency style plot:
plot(sort(x, decreasing=TRUE), 1:length(x), log="xy",
     xlab=expression(paste("Values, ", italic(x))),
     ylab=expression( paste("Number of ", values >= x)), mgp=mgpVals)
  # , xlim=c(2, 500), ylim=c(0.8, 40000), axes=FALSE, xaxs="i", yaxs="i", mgp=mgpVals)

x.PLB = seq(min(x), max(x), length=1000)     # x values to plot PLB
y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                  xmax = max(x.PLB))) * length(x)    
lines(x.PLB, y.PLB, col="red") #, lty=5)

legend("topright", paste("(h) MLE exponent=", signif(PLB.bMLE, 3)), bty="n",
       inset=inset)

# To add the curves at the limits of the 95% confidence interval:
#for(i in c(1, length(bIn95)))       # for(i in 1:length(bIn95))  to see all vals
#    {
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#                  xmax = max(x.PLB))) * length(x), col="red", lty=3)
#    }  

dev.off()


# Standard histogram, but with a break in the y-axis. See PLBfunctions.r for
#  code.
postscript("standHist10000.eps", height = 2.7, width = figwidth/2,
           horizontal=FALSE,  paper="special")  #  height=4, width=4   was 7,4
par(mai=c(0.6, 0.6, 0.15, 0.3))
# mgpVals = c(2, 0.5, 0)            # mgp values   2.0, 0.5, 0
gap.barplot.cust(hLlin.list$counts, gap=c(9,988),
                 ytics = c(seq(0, 8, by=4), seq(988, 996, by=4)),
                 midpoints = hLlin.list$mids,
                 breakpoints = hLlin.list$breaks,
                 xlim = c(0,100*ceiling(xmax/100)),
                 yaxs="i",
                 ylim= c(0, 17),  # max = max(y) - gap[2] + gap[1] + some space
                 col=rep("black", length(hLlin.list$counts)),
                 xaxs="i", xlab=expression(paste("Values, ", italic(x))),
                 ylab="Count in each bin", mgp=c(1.8,0.5,0))
                               # , default: mgp=c(3,1,0))# , mgp=c(1.8,0.5,0))
# mgp is margin line in mex units, for axis title, axis labels and axis line
dev.off()

save.image(file = "fitting2-10000.RData")


