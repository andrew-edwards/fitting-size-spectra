# fitting3rep-bMinus25.r - redoing with b=-2.5. 2/12/15.

# fitting3rep-bMinus15.r - redoing with b=-1.5. 2/12/15.

# fitting3rep.r - as for fitting2rep.r, but now plotting histograms in terms
#  of estimates of b rather than estimates of slope, since they will be clearer.
#  13th July 2015.

# fitting2rep.r - repeatedly fitting simulated data sets (drawn from PLB
#  with known fixed parameters) using the various techniques, and plotting
#  histograms of estimated slopes/exponents. Based on fitting2.r and
#  some of fitting1repr. 3rd July 2015

# fitting1rep.r - repeatedly fitting simulated data sets (drawn from PLB
#  with known fixed parameters) using the various techniques, and plotting
#  histograms of estimated slopes/exponents. 10th September 2014.

rm(list=ls())
print(date())
# source("PLBfunctions.r")  # to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy). Must come after load(..RData)
                          #  else functions can get overwritten.
redo.simulation = FALSE   # Whether or not to redo the simulations, 0 or 1
if(!redo.simulation)
  {load("fitting3rep-bMinus25.RData")
  source("../PLBfunctions.r")
  } else
  {                                                   
source("../PLBfunctions.r")
n = 1000                  # sample size
b.known = -2.5              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 1000         # known fixed value of xmax

num.reps = 10000          # number of times to draw sets of n random numbers.
                          #  (throwing n PLB dice num.reps times)  
set.seed(42)              # Same seed as for fitting1.r. 

# Record the slope for each method
Llin.rep = numeric(num.reps)*NA
LT.rep = Llin.rep
LTplus1.rep = Llin.rep
LBmiz.rep = Llin.rep 
LBbiom.rep = Llin.rep
LBNbiom.rep = Llin.rep
LCD.rep = Llin.rep 
MLE.rep = Llin.rep
MLEfix.rep = Llin.rep  # Adding in MLE calculations when we fix xmax=xmax.known

# Record the confidence intervals (in hindsight could have maybe done lists)
Llin.rep.conf = data.frame(confMin = Llin.rep, confMax = Llin.rep)
LT.rep.conf = Llin.rep.conf
LTplus1.rep.conf = Llin.rep.conf
LBmiz.rep.conf = Llin.rep.conf
LBbiom.rep.conf = Llin.rep.conf
LBNbiom.rep.conf = Llin.rep.conf
LCD.rep.conf = Llin.rep.conf 
MLE.rep.conf = Llin.rep.conf
MLEfix.rep.conf = Llin.rep.conf

MLE.rep.xmax = Llin.rep   # Also save the xmax =max(x) for each run, to see how
                           #  correlates with estimate of b.

num.bins = 8    # number of bins for standard histogram and Llin method, though
                #  this is only a suggestion (and can get overridden). Daan used
                #  8 bins.
hLBmiz.num.bins = num.bins    # for mizer method
# Main loop for doing the fitting num.reps times
for(iii in 1:num.reps)
{
if(iii %in% seq(1000, num.reps, 1000)) paste("iii = ", iii)    # show progress

x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)

log.x = log(x)                      # to avoid keep calculating
sum.log.x = sum( log.x ) 
xmin = min(x)
xmax = max(x)

# Notation:
# hAAA - h(istrogram) for method AAA.

# Llin method - plotting binned data on log-linear axes then fitting regression,
#  as done by Daan et al. 2005.
hLlin.list = Llin.method(x, num.bins = num.bins)
Llin.rep[iii] = hLlin.list$slope
Llin.rep.conf[iii,] = hLlin.list$confVals


# LT method - plotting binned data on log-log axes then fitting regression,
#  as done by Boldt et al. 2005, natural log of counts plotted against natural
#  log of size-class midpoints.

# Use Llin method's binning.
hLT.log.mids = log(hLlin.list$mids)
hLT.log.counts = log(hLlin.list$counts)
hLT.log.counts[ is.infinite(hLT.log.counts) ] = NA  
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

hLT.lm = lm( hLT.log.counts ~ hLT.log.mids, na.action=na.omit)
LT.rep[iii] = hLT.lm$coeff[2]
LT.rep.conf[iii,] = confint(hLT.lm, "hLT.log.mids", 0.95)

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
LTplus1.rep[iii] = hLTplus1.lm$coeff[2]
LTplus1.rep.conf[iii,] = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95)

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

beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est

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
LBmiz.rep[iii] = hLBmiz.lm$coeff[2]
LBmiz.rep.conf[iii,] = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95)

# LBbiom method - binning data using log2 bins, calculating biomass (not counts)
#  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
#  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

hLBNbiom.list = LBNbiom.method(x)    # Does this method and the next.
LBbiom.rep[iii] = hLBNbiom.list[["unNorm.slope"]]
LBbiom.rep.conf[iii, ] = hLBNbiom.list[["unNorm.conf"]]

# LBNbiom method - on biomass, not counts, as per Julia's 2005 paper.
#  log2 bins of bodymass, sum the total biomass in each bin, normalise
#  biomasses by binwidths, fit regression to log10(normalised biomass) v
#  log10(midpoint of bin).

# hLBNbiom.list = LBNbiom.method(x) - already done above

LBNbiom.rep[iii] = hLBNbiom.list[["norm.slope"]]
LBNbiom.rep.conf[iii, ] = hLBNbiom.list[["norm.conf"]]

# Cumulative Distribution, LCD method
x.sorted = sort(x, decreasing=TRUE)
logSorted = log(x.sorted)
logProp = log((1:length(x))/length(x))

hLCD.lm = lm(logProp ~ logSorted)   # plot(fitsortedlog10) shows
                                                 #  residuals not good
LCD.rep[iii] = hLCD.lm$coeff[2]
LCD.rep.conf[iii,] = confint(hLCD.lm, "logSorted", 0.95)

# MLE (maximum likelihood method) calculations. Estimate xmax=max(x)

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
#  as a starting point for nlm for MLE of b for PLB model.
PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1
    
PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
    xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )

PLB.bMLE = PLB.minLL$estimate
MLE.rep[iii] = PLB.bMLE
MLE.rep.xmax[iii] = xmax

# 95% confidence intervals for MLE method.
PLB.minNegLL = PLB.minLL$minimum

# Values of b to test to obtain confidence interval. For the movement data
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
MLE.rep.conf[iii,] = c(PLB.MLE.bConf[1], PLB.MLE.bConf[2])

# MLE (maximum likelihood method) calculations, but fix xmax=xmax.known

PLBfix.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
    xmin=xmin, xmax=xmax.known, sumlogx=sum.log.x) #, print.level=2 )

PLBfix.bMLE = PLBfix.minLL$estimate
MLEfix.rep[iii] = PLBfix.bMLE

# 95% confidence intervals for MLE method.
PLBfix.minNegLL = PLBfix.minLL$minimum

# Values of b to test to obtain confidence interval. For the movement data
#  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
#  symmetric interval here.
# bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.001) 
    
PLBfix.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec))
    {
        PLBfix.LLvals[i] = negLL.PLB(bvec[i], x=x, n=length(x), xmin=xmin,
            xmax=xmax.known, sumlogx=sum.log.x)   
    }
critVal = PLBfix.minNegLL  + qchisq(0.95,1)/2
                    # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLBfix.LLvals < critVal ]
                    # b values in 95% confidence interval
PLBfix.MLE.bConf = c(min(bIn95), max(bIn95))
if(PLBfix.MLE.bConf[1] == min(bvec) | PLBfix.MLE.bConf[2] == max(bvec))
  { windows()
    plot(bvec, PLBfix.LLvals)
    abline(h = critVal, col="red")
    stop("Need to make bvec larger for PLBfix - see R window")   # Could automate
  }
MLEfix.rep.conf[iii,] = c(PLBfix.MLE.bConf[1], PLBfix.MLE.bConf[2])


}  # End for for(iii in 1:num.reps) loop

} # End of if(!redo.simulation) {load("fitting1rep.RData")} else { 


# brange = c(-3.5,0) #c(-2.25, 0) # -0.75)    # Histograms of results
# brange = range( c( Llin.rep, LT.rep, LTplus1.rep, LBmiz.rep,
#    LBbiom.rep, LBNbiom.rep, MLE.rep))
# Comes out as -3.38, 0.233 for seed=42, num.reps = 10000
xrange = c(-3.5, 0.5)-0.5         # range of x-axis for histograms - actually to define bins
# xbigticks = -c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25)    # CHANGE TO DEPEND ON brange
xbigticks = seq(-3.5, 0.5, 0.5)

xsmallticks = seq(xrange[1], xrange[2], by=0.1)
# breakshist = seq(xrange[1]-0.05/2, xrange[2]+0.05, by=0.05)   # then 2 is centre of a bin, but then it missed a value 1.001   (what does that mean?)
# breakshist = seq(xrange[1], xrange[2], by=0.05)
         # ends a bin at 0 for xrange[1]=0, just won't have 2 in centre of bin
breakshist = seq(xrange[1], xrange[2], by=3/61)
         # bin ends at 0, 2 is centered. width of 4/81 from solving
         #  the mean of a bin equals 2, ((N+1)w + Nw)/2 = 2, and setting number
         #  of the bin that includes 2, N, =40 (which is how many you get with
         #  0.05, which looked good).

# xrange=c(-3.5, 0.5), think I need:
#   -2 to be a midpoint, which is 1.5 above the minimum. 30th bin contains -2,
#   so solve ((N+1)w + Nw)/2 = 1.5 with N = 30, gives 61w/2 = 3/2, w = 3/61

# breakshist = seq(xrange[1], xrange[2], by=0.025)
figheight = 7 # 5.6
figwidth = 5.7    # 5.7 inches for JAE
postscript("fitting3rep-bMinus25.eps", height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")  
par(omi = c(0.12, 0.05, 0.12, 0.0))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size 
par(mai=c(0.5, 0.5, 0.1, 0.3))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

vertCol = "red"            # Colour for vertical lines
vertThick = 1              # Thickness for vertical lines

# ylimA = c(0, num.reps)     # Want to define just one y range if possible.
xLim = c(-3.5, -1.8)    # c(-3.2, -1.5)
xLimLlin = c(-1.2, 0.5) # same as b=-2
ylimA = c(0, 8000)
ylimLlin = c(0, 10000)
# ylimB = c(0, 3250)
cexAxis = 0.9      # font size for axes labels to make the y ones fit okay

# Llin.rep has different breakhist so that 0 is a breakpoint. Same width as
#  others, just shifted. 
hist(Llin.rep, xlim=xLimLlin, breaks=breakshist - breakshist[31],
  xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimLlin)  #  ylim=c(0,1100), 
axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks
xsmallticksLlin = seq(xLimLlin[1], xLimLlin[2], by=0.1)
axis(1, at=xsmallticksLlin, labels=rep("",length(xsmallticksLlin)), tcl=-0.2)
axis(2, at=c(0, 5000, 10000),
         labels = c(0, 5000, 10000),
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
axis(2, at=seq(0, 10000, 1000),
         labels = rep("", 11),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
# axis(2, at=seq(0, 10000, 500),
#          labels = rep("", 14), mgp=c(1.7,0.7,0), tcl=-0.2)  # small ticks
abline(v=b.known, col=vertCol, lwd=vertThick)
inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(a) Llin", bty="n", inset=inset)
#axis(2, at=seq(0, 6500, 500),
#        labels = c(0, "", "", "", 2000, "", "", "", 4000, "", "", "", 6000, ""),
#         mgp=c(1.7,0.7,0))  # big ticks
# axis(2, at=seq(0, 1100, by=100), labels = rep("", 12), mgp=c(1.7,0.7,0), tcl=-0.2)  # small ticks

# replace with legend
figlabpos = 0.93    # proportion of x and y axis lengths to put (a) in etc.
# xlabpos = 1.5 * (1 - figlabpos) *0.7 + 0.75 # * 0.8 as axis so stretched,=0.82
xlabpos = 0.75       # just play with a number, as now using pos=4 in text
text( xlabpos, figlabpos * 1100, "(a) Llin", pos=4)


hist(LT.rep, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(b) LT", bty="n", inset=inset)

hist(LTplus1.rep, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(c) LTplus1", bty="n", inset=inset)

hist(LBmiz.rep-1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(d) LBmiz", bty="n", inset=inset)

hist(LBbiom.rep-2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(e) LBbiom", bty="n", inset=inset)

hist(LBNbiom.rep-1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(f) LBNbiom", bty="n", inset=inset)


hist(LCD.rep-1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(g) LCD", bty="n", inset=inset)

hist(MLE.rep, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(h) MLE", bty="n", inset=inset)

# text(-2.5, 0, expression(paste("Estimate of ", italic(b))))
mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)
# mtext("hello", side=1, outer=TRUE, line=-1)

dev.off()

# Histograms for MLE and MLEfix methods.
postscript("fitting3repMLEfix-bMinus25.eps", height = 0.8*figheight, width = 0.8*figwidth,
           horizontal=FALSE,  paper="special")  
par(omi = c(0.12, 0.05, 0.12, 0.0))      # outer margins in inches
par(mfrow=c(2,1)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size 
par(mai=c(0.5, 0.5, 0.1, 0.3))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

# xrangeMLEfix = c(-2.225, -1.775)  # Range for these two figures
# breakshistMLEfix = seq(xrangeMLEfix[1], xrangeMLEfix[2], by=3/61)

xLimFix = c(-2.5, -1.5)
yLimFix = c(0, 6000)
hist(MLE.rep, xlim=xLimFix, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = yLimFix)  #  ylim=c(0,1100),
histAxes()
axis(2, at=6000, labels = 6000, mgp=c(1.7,0.7,0), cex.axis=cexAxis)
                                        # big tick label
legend("topleft", "(a) MLE", bty="n", inset=0.5*inset)

hist(MLEfix.rep, xlim=xLimFix, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = yLimFix)  #  ylim=c(0,1100),
histAxes()
axis(2, at=6000, labels = 6000, mgp=c(1.7,0.7,0), cex.axis=cexAxis)
                                        # big tick label
legend("topleft", "(b) MLEfix", bty="n", inset=0.5*inset)

# text(-2.5, 0, expression(paste("Estimate of ", italic(b))))
mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)
# mtext("hello", side=1, outer=TRUE, line=-1)

dev.off()



# Prints Latex code for table that summarises the results
# Going to do in terms of b, not slope, so subtract 1 or 2 from some results:
# Quartiles:
print("Method & Slope is estimating: & 1st quartile & Median & Mean & 3rd quartile & Percentage below true \\")
print(paste("Llin & - & ", qqtab(Llin.rep, dig=3), "\\"))
print(paste("LT & $b$ & ", qqtab(LT.rep), "\\"))
print(paste("LTplus1 & $b$ & ", qqtab(LTplus1.rep), "\\"))
print(paste("LBmiz & $b+1$ & ", qqtab(LBmiz.rep - 1), "\\"))
print(paste("LBbiom & $b+2$ & ", qqtab(LBbiom.rep - 2), "\\"))
print(paste("LBNbiom & $b+1$ & ", qqtab(LBNbiom.rep - 1), "\\"))
print(paste("LCD & $b+1$ & ", qqtab(LCD.rep - 1), "\\"))
print(paste("MLE & $b$ & ", qqtab(MLE.rep), "\\hline"))
# To give slopes, rather than converting them to b:
#print(paste("LBmiz & ", qqtab(LBmiz.rep, true=b.known+1), "\\"))
#print(paste("LBbiom & ", qqtab(LBbiom.rep, true=b.known+2), "\\"))
#print(paste("LBNbiom & ", qqtab(LBNbiom.rep, true=b.known+1), "\\"))
#print(paste("LCD & ", qqtab(LCD.rep, true=b.known+1), "\\"))

# 5% and 95% values:
print("Method & Slope represents & 5% quantile & Median & Mean & 95% quantile & Percentage below true \\")
print(paste("Llin & - & ", qqtab(Llin.rep, dig=2, quants=c(0.05, 0.95)), "\\"))
print(paste("LT & $b$ & ", qqtab(LT.rep, quants=c(0.05, 0.95)), "\\"))
print(paste("LTplus1 & $b$ & ", qqtab(LTplus1.rep, quants=c(0.05, 0.95)), "\\"))
print(paste("LBmiz & $b+1$ & ", qqtab(LBmiz.rep - 1, quants=c(0.05, 0.95)), "\\"))
print(paste("LBbiom & $b+2$ & ", qqtab(LBbiom.rep - 2, quants=c(0.05, 0.95)), "\\"))
print(paste("LBNbiom & $b+1$ & ", qqtab(LBNbiom.rep - 1, quants=c(0.05, 0.95)), "\\"))
print(paste("LCD & $b+1$ & ", qqtab(LCD.rep - 1, quants=c(0.05, 0.95)), "\\"))
print(paste("MLE & $b$ & ", qqtab(MLE.rep, quants=c(0.05, 0.95)), "\\hline"))
print(paste("% MLEfix & $b$ & ", qqtab(MLEfix.rep, quants=c(0.05, 0.95)), "\\"))
# For slopes:
#print(paste("LBmiz &  ", qqtab(LBmiz.rep, true=b.known+1, quants=c(0.05, 0.95)), "\\"))
#print(paste("LBbiom &  ", qqtab(LBbiom.rep, true=b.known+2, quants=c(0.05, 0.95)), "\\"))
#print(paste("LBNbiom & ", qqtab(LBNbiom.rep, true=b.known+1, quants=c(0.05, 0.95)), "\\"))
#print(paste("LCD & ", qqtab(LCD.rep, true=b.known+1, quants=c(0.05, 0.95)), "\\"))


save.image(file = "fitting3rep-bMinus25.RData")
