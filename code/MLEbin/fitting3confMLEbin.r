# fitting3confMLEbin.r - histogram of MLE values and confidence intervals for
#  the repeated simulated binned data sets, using the MLEbin method. To be new
#  figure in revised version of first manuscript.
#  3rd May 2015.

# fitting3confBin.r - histogram of MLE values and confidence intervals for
#  the repeated simulated binned data sets, not taking into account the binning
#  (MLEmid method) and taking binning into account (MLEbin method).
#  22th October 2015 (the future).

# fitting3conf.r - confidence intervals based on b, not slopes. 13th July 2015

# fittingConf.r - showing the confidence intervals for some of
#  methods, and how many of them include the true value. 8th July 2015

rm(list=ls())
require(dplyr)

load("fitting3repMLEbin.RData")
source("../PLBfunctions.r")  # to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy)

figheight = 4.1 # 4x2 are 7 
figwidth = 5.7/2    # 5.7 inches for JAE

inCol = "black"           # Colour for the true value being within the 95% CI
outCol = "red"            # Colour for the true value being outside the 95% CI

postscript("fitting3confMLEbin.eps", height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")  

par(omi = c(0.14, 0, 0.1, 0.20))      # outer margins in inches
par(mfrow=c(2,1)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size 
par(mai=c(0.3, 0.5, 0.08, 0))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small
vertThick = 1              # Thickness for vertical lines

# Histogram of MLEs

vertCol = "red"            # Colour for vertical lines
vertThick = 1              # Thickness for vertical lines

# range(MLE.rep.conf)      -2.078214 -1.753170
# range(MLEbin.rep.conf)   -2.209575 -1.823068
# xrange = c(-2.25, -1.75)     # common width of axes, symmetric about -2

# from fitting3rep.r, make x-range the same
xrange = c(-3.5, 0.5)         # range of x-axis for histograms 
xbigticks = seq(-3.5, 0.5, 0.5)

xsmallticks = seq(xrange[1], xrange[2], by=0.1)
# breakshist = seq(xrange[1], xrange[2], by=0.05)
         # ends a bin at 0 for xrange[1]=0, just won't have 2 in centre of bin
breakshist = seq(xrange[1], xrange[2], by=3/61)

xLim = c(-3.2, -1.5)
ylimA = c(0, 5500)

cexAxis = 1  # was 0.9      # font size for axes labels to make y ones fit ok

inset = c(-0.08, -0.04)       # (default in confPlot)
xlabpos = 0.75       # just play with a number, as now using pos=4 in text

hist(MLEbin.rep, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100), 
histAxes()
# Can prob delete:
# It only did 2000, so add extra ones manually:
# axis(2, at=c(0, 1000, 3000),    
#          labels = c(0, 1000, 3000),     
#         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled

legend("topleft", "(a) MLEbin", bty="n", inset=inset)

# Each of these plots a panel for one method. Define xLim if the default
#  (integer-based calculation) is not suitable

#MLEmid.rep.conf.sort = confPlot(MLE.rep.conf, legName="(c) MLEmid",
#    xLim = xrange, xsmallticks = xsmallticks, insetVal=inset,
#    insetVal2=inset + c(0, 0.1), yLabels=FALSE)
#axis(2, at = c(0, 100, 200, 300), tck=-0.04)  # confPlot did 50, 150,..

xLimCom = c(-2.65, -1.5)     # common width of axes for most methods in main
                             #  figure

MLEbin.rep.conf.sort = confPlot(MLEbin.rep.conf, legName="(b) MLEbin",
    xLim = xLimCom, xsmallticks = xsmallticks, insetVal=inset,
    insetVal2=inset + c(0, 0.1), yLabels = FALSE)
axis(2, at = c(0, 100, 200, 300), tck=-0.04)  # confPlot did 50, 150,.. 

mtext(expression(paste("Estimate of ", italic(b)), sep=""),
       side=1, outer=TRUE, line=-0.2, cex=0.8)

dev.off()

print("Range of widths of confidence intervals for MLEbin method is")
print(range(MLEbin.rep.conf[2] - MLEbin.rep.conf[1]))
