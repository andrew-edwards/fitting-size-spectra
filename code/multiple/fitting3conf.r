# fitting3conf.r - confidence intervals based on b, not slopes. 13th July 2015

# fittingConf.r - showing the confidence intervals for some of
#  methods, and how many of them include the true value. 8th July 2015

# fitting2rep.r - repeatedly fitting simulated data sets (drawn from PLB
#  with known fixed parameters) using the various techniques, and plotting
#  histograms of estimated slopes/exponents. Based on fitting2.r and
#  some of fitting1repr. 3rd July 2015

# fitting1rep.r - repeatedly fitting simulated data sets (drawn from PLB
#  with known fixed parameters) using the various techniques, and plotting
#  histograms of estimated slopes/exponents. 10th September 2014.


rm(list=ls())
require(dplyr)

load("fitting3rep.RData")
source("../PLBfunctions.r")# to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy)

inCol = "black"           # Colour for the true value being within the 95% CI
outCol = "red"            # Colour for the true value being outside the 95% CI

postscript("fitting3conf.eps", height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")  

par(omi = c(0.14, 0, 0.1, 0.15))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size 
par(mai=c(0.3, 0.5, 0.08, 0))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small
vertThick = 1              # Thickness for vertical lines

# Each of these plots a panel for one method. Define xLim if the default
#  (integer-based calculation) is not suitable
Llin.rep.conf.sort = confPlot(Llin.rep.conf, legName="(a) Llin",
    xLim = c(-0.25, 0.25)) #, colourCode=FALSE)

LT.rep.conf.sort = confPlot(LT.rep.conf, legName="(b) LT", yLab="", yLabels=FALSE)

LTplus1.rep.conf.sort = confPlot(LTplus1.rep.conf, legName="(c) LTplus1")

xLimCom = c(-2.65, -1.5)     # common width of axes for the remainder

# range(LBmiz.rep.conf) =  -1.5894281 -0.5268428. Was -1.6, -0.4
LBmiz.rep.conf.sort = confPlot(LBmiz.rep.conf-1, legName="(d) LBmiz",
    xLim = xLimCom, yLab="", yLabels=FALSE)

# range(LBbiom.rep.conf) = -0.6106510  0.4370617
LBbiom.rep.conf.sort = confPlot(LBbiom.rep.conf-2, legName="(e) LBbiom",
      xLim = xLimCom)

# range(LBNbiom.rep.conf) = -1.6106510 -0.5629383
LBNbiom.rep.conf.sort = confPlot(LBNbiom.rep.conf-1, legName="(f) LBNbiom",
     yLab="", yLabels=FALSE, xLim = xLimCom)

# range(LCD.rep.conf) = -1.1744680 -0.8713301   was c(-1.2, -0.8)
LCD.rep.conf.sort = confPlot(LCD.rep.conf-1, legName="(g) LCD",
    xLim = xLimCom, vertFirst = TRUE)

# c(-2.2, -1.8) - looks good, but better to have consistent
MLE.rep.conf.sort = confPlot(MLE.rep.conf, legName="(h) MLE",
    xLim = xLimCom, yLab="", yLabels=FALSE)
mtext(expression(paste("Estimate of ", italic(b)), sep=""),
      side=1, outer=TRUE, line=-0.2, cex=0.8)

dev.off()

# Confidence intervals for MLE and MLEfix methods.
postscript("fitting3confMLEfix.eps", height = 0.8*figheight,
           width = 0.8*figwidth,
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

MLE.rep.conf.sort = confPlot(MLE.rep.conf, legName="(a) MLE",
    xLim = xLimCom, insetVal = c(-0.04, -0.03), insetVal2 = c(-0.04, 0.05))
                                        

MLEfix.rep.conf.sort = confPlot(MLEfix.rep.conf, legName="(b) MLEfix",
    xLim = xLimCom, insetVal = c(-0.04, -0.03), insetVal2 = c(-0.04, 0.05))

mtext(expression(paste("Estimate of ", italic(b)), sep=""),
      side=1, outer=TRUE, line=-0.2, cex=0.8)

dev.off()

print("Range of widths of confidence intervals for MLE method is")
print(range(MLE.rep.conf[2] - MLE.rep.conf[1]))

print("Range of widths of confidence intervals for MLEfix method is")
print(range(MLEfix.rep.conf[2] - MLEfix.rep.conf[1]))
