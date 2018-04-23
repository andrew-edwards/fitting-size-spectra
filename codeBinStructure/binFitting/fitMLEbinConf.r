# fitMLEbinConf.r - histogram of MLE values and confidence intervals for
#  the repeated simulated binned data sets, using the MLEmid and MLEbin
#  method(s) for four binning types. Plots 2x2 figure for each binning
#  type, and 4x2 figures of MLEs and of confidence intervals for all
#  binning types and both methods.
#  Adapting from fitting3confMLEbin.r from first manuscript.
#  6th June 2016.

# fitting3confMLEbin.r - histogram of MLE values and confidence intervals for
#  the repeated simulated binned data sets, using the MLEbin method. To be new
#  figure in revised version of first manuscript.
#  3rd May 2015.

# fitting3confBin.r - histogram of MLE values and confidence intervals for
#  the repeated simulated binned data sets, not taking into account the binning
#  (MLEmid method) and taking binning into account (MLEbin method).
#  22th October 2015 (the future).


rm(list=ls())
require(dplyr)

load("fitMLEmidMLEbin.RData")
source("../../code/PLBfunctions.r")
source("../countsFunctions.r")
                          # to reload in required functions in case changed
if(MLEmethods > 2) stop("Need to adapt code to deal with >1 MLE method")

figheight = 4.1 # 4x2 are 7
figwidth = 5.7    # 5.7 inches for JAE

vertCol = "red"            # Colour for vertical lines in MLE histogram
vertThick = 1              # Thickness for vertical lines

inCol = "black"           # Colour for the true value being within the 95% CI
outCol = "red"            # Colour for the true value being outside the 95% CI

print("Overall range of MLE values:")
print(range(MLE.array[ , , ]))
print("Overall range of confidence intervals:")
print(range(MLEconf.array[ , , , ]))

# Common x axis and bin breaks for histograms
xrange = c(-2.4, -1.1)     # common width of axes
                            # Need to play around
                            #  with but do for now
xbigticks = seq(xrange[1], xrange[2], by=0.4)
xsmallticks = seq(xrange[1], xrange[2], by=0.1)

yBigTickLab = seq(0, num.reps, 3000)
yBigTickNoLab = seq(0, num.reps, 1000)
ySmallTick = seq(0, num.reps, 500)

# Want -2 to be a midpoint of the Nth bin, which is yy above the minimum.
#  Say the 20th bin contains -2, so solve ((N+1)w + Nw)/2 = yy
#  gives  w = 2 * yy / (2*N + 1)
binwidth = 2 * (b.known - xrange[1]) / ( 2 * 10 + 1)
breakshist =  seq(xrange[1], length=ceiling((xrange[2] - xrange[1])/binwidth)+1,
    by=binwidth)

ylimA = c(0, 7000)

for(binTypeInd in 1:binTypes)
  {
  postscript(paste0("fitMLEbinConf", binType[[binTypeInd]], ".eps"),
           height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")

  par(omi = c(0.14, 0, 0.1, 0.15))      # outer margins (inches; set last to 0.2
                                        #  for single column)
  par(mfrow=c(2,2))

  oldmai = par("mai")
  par(mai=c(0.3, 0.5, 0.08, 0))  # Affects all four figures if don't change again
  par(xaxs="i", yaxs="i")    # Have to define here for hist
  par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
  par(cex = 0.8)             # With no option all text comes out a bit small

  # Histogram of MLEs for the two methods

  cexAxis = 1                # font size for axes labels to make the y ones fit ok

  inset = c(0, -0.04)        # (default in confPlot)
  # xlabpos = 0.75       # just play with a number, as now using pos=4 in text


  hist(MLE.array[ ,binTypeInd, "MLEmid"], xlim=xrange, breaks=breakshist,
       xlab="", ylab="Frequency", main="", axes=FALSE, ylim = ylimA)
  histAxes(yBigTickLab = yBigTickLab, yBigTickNoLab = yBigTickNoLab,
           ySmallTick = ySmallTick)
  legend("topright", "(a) MLEmid", bty="n", inset=inset)

  hist(MLE.array[ ,binTypeInd, "MLEbin"], xlim=xrange, breaks=breakshist,
       xlab="", ylab="Frequency", main="", axes=FALSE, ylim = ylimA)
  histAxes(yBigTickLab = yBigTickLab, yBigTickNoLab = yBigTickNoLab,
           ySmallTick = ySmallTick)
  legend("topright", "(b) MLEbin", bty="n", inset=inset)


  # Each of these plots a panel for one method. Define xLim if the default
  #  (integer-based calculation) is not suitable

  MLEmid.rep.conf.sort = confPlot(as.data.frame(
    MLEconf.array[ ,binTypeInd, "MLEmid", ]),
    legName="(c) MLEmid",
    xLim = xrange, xsmallticks = xsmallticks, insetVal=inset,
    insetVal2=inset + c(0, 0.1), yLabels=FALSE, legLoc="topright")

  axis(2, at = c(0, 100, 200, 300), tck=-0.04)  # confPlot did 50, 150,..

  MLEbin.rep.conf.sort = confPlot(as.data.frame(
      MLEconf.array[ ,binTypeInd, "MLEbin", ]),
      legName="(d) MLEbin",
      xLim = xrange, xsmallticks = xsmallticks, insetVal=inset,
      insetVal2=inset + c(0, 0.1), yLabels = FALSE, legLoc="topright")
  axis(2, at = c(0, 100, 200, 300), tck=-0.04)  # confPlot did 50, 150,..

  mtext(expression(paste("Estimate of ", italic(b)), sep=""),
       side=1, outer=TRUE, line=-0.2, cex=0.8)

  dev.off()
}                                   # End for(binTypeInd in 1:binTypes)



figheight = 7     # 4x2 figs are 7 inches
figwidth = 5.7

# Histograms of the MLE values for both methods and all four binning types.
postscript("fitMLEhists.eps",
           height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")

par(omi = c(0.12, 0.05, 0.2, 0.0))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

par(mai=c(0.5, 0.5, 0.0, 0.3))
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

cexAxis = 0.9              # font size for axes labels to make the y ones
                           #  fit okay
inset = c(0, -0.04)

legLabMid = c("(a)", "(c)", "(e)", "(g)")    # label each panel in MLEmid column
legLabBin = c("(b)", "(d)", "(f)", "(h)")    # label each panel in MLEmid column

for(binTypeInd in 1:binTypes)
  {
  hist(MLE.array[ ,binTypeInd, "MLEmid"], xlim=xrange, breaks=breakshist,
         xlab="", ylab="", main="", axes=FALSE, ylim = ylimA)
  histAxes(yBigTickLab = yBigTickLab, yBigTickNoLab = yBigTickNoLab,
           ySmallTick = ySmallTick)
  legend("topright", paste(legLabMid[binTypeInd], binType.name[binTypeInd]),
                           bty="n", inset=inset)

  hist(MLE.array[ ,binTypeInd, "MLEbin"], xlim=xrange, breaks=breakshist,
       xlab="", ylab="", main="", axes=FALSE, ylim = ylimA)
  histAxes(yBigTickLab = yBigTickLab, yBigTickNoLab = yBigTickNoLab,
           ySmallTick = ySmallTick)
  legend("topright", paste(legLabBin[binTypeInd], binType.name[binTypeInd]),
                           bty="n", inset=inset)
}                # End for(binTypeInd in 1:binTypes) for 4x2 MLE hist figure.

mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)
mtext("Frequency", side=2, outer=TRUE, line=-1)

mtext("    MLEmid                                                MLEbin",
      side=3, outer=TRUE, line=0)

dev.off()

# Confidence intervals for the MLEs for both methods and all four binning types.
postscript("fitMLEconfs.eps",
           height = figheight, width = figwidth,
           horizontal=FALSE,  paper="special")

par(omi = c(0.2, 0.05, 0.22, 0.1))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

par(mai=c(0.3, 0.5, 0.08, 0))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small
vertThick = 1              # Thickness for vertical lines

# Need a different inset for panel (e) MLEmid Linear 10.
insetMat = matrix(rep(c(-0.01, -0.04), 4), ncol=2, byrow=TRUE)
insetMat[3, 1] = 0.3

for(binTypeInd in 1:binTypes)
  {
  # confPlot returns a data.frame of intervals, but no need to save them
  #  for each binType.
  res = confPlot(as.data.frame(
    MLEconf.array[ ,binTypeInd, "MLEmid", ]),
    legName=paste(legLabMid[binTypeInd], binType.name[binTypeInd]),
    xLim = xrange, xsmallticks = xsmallticks, insetVal=insetMat[binTypeInd,],
    insetVal2=insetMat[binTypeInd,] + c(0, 0.12), legLoc="topright", yLab="")

  res = confPlot(as.data.frame(
    MLEconf.array[ ,binTypeInd, "MLEbin", ]),
    legName=paste(legLabBin[binTypeInd], binType.name[binTypeInd]),
    xLim = xrange, xsmallticks = xsmallticks, insetVal=inset,
    insetVal2=inset + c(0, 0.12), yLabels=FALSE, legLoc="topright", yLab="")
}                # End for(binTypeInd in 1:binTypes) for 4x2 conf figure.

mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=0)
mtext("Sample number", side=2, outer=TRUE, line=-1)

mtext("         MLEmid                                              MLEbin",
      side=3, outer=TRUE, line=0)

dev.off()

# 5% and 95% values of MLEs, to copy into .tex write up.
print("Binning type & Method & 5% quantile & Median & Mean & 95% quantile & Percentage below true \\ \\hline", quote=FALSE)

for(j in 1:binTypes)                    # Loop over binning type
  {
  print(paste(binType.name[[j]], " & ", MLEmethod.name[1], " & ",
      qqtab(MLE.array[ ,j,1], quants=c(0.05, 0.95)), "\\"), quote=FALSE)
  print(paste(" & ", MLEmethod.name[2], " & ",
      qqtab(MLE.array[ ,j,2], quants=c(0.05, 0.95)), "\\"), quote=FALSE)
}

# Range of widths of confidence intervals, may want to quote them but
#  unlikely need to report the full table.
print("Binning type & Method & Shortest CI & Widest CI\\ \\hline", quote=FALSE)
for(binTypeInd in 1:binTypes)
  {
  rangeVals = range(MLEconf.array[ ,binTypeInd, "MLEmid", "confMax"] -
    MLEconf.array[ ,binTypeInd, "MLEmid", "confMin"])
  print(paste(binType.name[[binTypeInd]], " & ", MLEmethod.name[1], " & ",
    rangeVals[1], " & ", rangeVals[2], "\\"), quote=FALSE)

  rangeVals = range(MLEconf.array[ ,binTypeInd, "MLEbin", "confMax"] -
    MLEconf.array[ ,binTypeInd, "MLEbin", "confMin"])
  print(paste(" & ", MLEmethod.name[2], " & ",
    rangeVals[1], " & ", rangeVals[2], "\\"), quote=FALSE)
}

