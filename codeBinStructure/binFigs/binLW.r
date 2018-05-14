# binLW.r - produce (automated) length-weight curves for two species.

require(dplyr)
# require(gplots)                 # for plotCI
require(shape)                  # for Arrows
rm(list=ls())

figheight = 5
figwidth = 5

rou = function(x, dig=1)
    {
    # Round a value to dig decimal places, and print that number of decimal
    #  places, even if it includes zero. round(79.01, 1) gives 79, but for
    #  Sweave we'd still often like it to be 79.0.
    # Args: x - quantity to be rounded
    #       dig - number of decimal places
    return(sprintf( paste("%.", dig, "f", sep=""),
                   round(x, dig=dig)))
    }

# Note that sp1 is Common Ling in the code, but then doing Lemon Sole first
#  in figure and text")

lengthToMass = function(lengths, LWa, LWb)
    {
        # Convert a vector of lengths to a vector of masses using the formula
        #  mass = LWa * lengths^LWb. LWa and LWb will be species-specific.
        # Args:
        #  lengths - vector of lengths
        #  LWa - multiplicative coefficient in equation
        #  LWb - exponent in equation
        if(min(c(lengths, LWa, LWb)) < 0)
            { stop("Need positive arguments in lengthToMass")  }
        return(LWa * lengths^LWb)
    }

sp1 = "Common Ling"   # Molva molva
sp2 = "Lemon Sole"    # Microstomus kitt
LWa = c(0.001, 0.0255)    # LIN then LEM   sp1 then sp2
LWb = c(3.4362, 2.7643)

col1 = c("red", "pink")           # Colours for species 1
col2 = c("blue", "lightblue")     # Colours for species 2
thick=7                           # thickness of bins

length = 10:50
length = 10:500                   # to check that the curves cross

mass = rbind(lengthToMass(length, LWa[1], LWb[1]),
             lengthToMass(length, LWa[2], LWb[2]))

# Length-weight relationships for two species, Figure 1.

postscript("binLW.eps", height = figheight, width = figwidth,
            horizontal=FALSE,  paper="special")

par(xaxs="i", yaxs="i")
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
lThick = 3                 # Thickness for curves
par(lend="butt")           # To have butted line caps, need for thick lines.

# Axes ranges not automated
plot(length, mass[1,], col=col1[1], type="l", xlab="Length, cm",
     ylab="Body mass, g", xlim=c(-6, 50), ylim=c(-60, 800), lwd=lThick)

lines(length, mass[2,], col=col2[1], lwd=lThick)
# Plotting up to length 200 and mass 50000 shows that the curves cross.

# Switching sp1 and sp2 around in the legend (and thus the text), since
#  looks more consistent.
legend("topleft", legend=c(sp2, sp1), lty=1, lwd=lThick, col=c(col2[1], col1[1]),
       bty="n", inset=c(0.1,-0.02))
axis(1, at = seq(0, 50, by=5), labels=rep("", 11), tck=-0.02)
axis(2, at = seq(0, 800, by=100), labels=rep("", 9), tck=-0.02)
axis(2, at = seq(0, 800, by=50), labels=rep("", 17), tck=-0.01)

lenBins = seq(10, 40, by=5)      # Length bin endpoints for the example
numBinBreaks = length(lenBins)
lenBinsStart = lenBins[-numBinBreaks]     # starting point for length bins
lenBinsEnd = lenBins[-1]

massBins = rbind(lengthToMass(lenBins, LWa[1], LWb[1]),
             lengthToMass(lenBins, LWa[2], LWb[2]))  # resulting species-specific
                                                     #  bins of body mass.
massBinsStart = massBins[, -numBinBreaks]
massBinsEnd = massBins[, -1]

colBins = rbind( rep(col1, length=round(length(lenBinsStart))),
    rep(col2, length=round(length(lenBinsStart))))
                                        # colours for bins for each species

yLengths = c(-40, -20)                  # y values to plot length bins, for
                                        #  each species
xMasses  = c(-3, -1)                    # x values to plot mass bins

# length bins:
segments(x0=lenBinsStart, y0=yLengths[1], x1 = lenBinsEnd, col=colBins[1,],
          lwd=thick)  # y1=y0
segments(x0=lenBinsStart, y0=yLengths[2], x1 = lenBinsEnd, col=colBins[2,],
          lwd=thick)

# resulting mass bins:
segments(x0=xMasses[1], y0=lengthToMass(lenBinsStart, LWa[1], LWb[1]),
         y1=lengthToMass(lenBinsEnd, LWa[1], LWb[1]),
         col=colBins[1,], lwd=thick)
segments(x0=xMasses[2], y0=lengthToMass(lenBinsStart, LWa[2], LWb[2]),
         y1=lengthToMass(lenBinsEnd, LWa[2], LWb[2]),
         col=colBins[2,], lwd=thick)


egBin = 5         # example bin number to highlight
offset = 0.1      # offset to shift species vertical lines so they show up
midOfBin = mean(lenBinsStart[c(egBin, egBin+1)])    # midpoint of example bin
massMidOfBin = rbind(lengthToMass(midOfBin, LWa[1], LWb[1]),
    lengthToMass(midOfBin, LWa[2], LWb[2]))

# example bin for species 1:
lines(c(rep(lenBinsStart[egBin], 2), xMasses[1])-offset,
       c(yLengths[1], rep(lengthToMass(lenBinsStart[egBin], LWa[1], LWb[1]), 2)),
       lty=3, lwd=1, col=col1[1])
lines(c(rep(lenBinsStart[egBin+1], 2), xMasses[1])-offset,
       c(yLengths[1], rep(lengthToMass(lenBinsStart[egBin+1], LWa[1], LWb[1]),
         2)), lty=3, lwd=1, col=col1[1])
lines(c(rep(midOfBin, 2), xMasses[1]) - offset,
       c(yLengths[1], rep(massMidOfBin[1,], 2)),
       lty=1, lwd=1, col=col1[1])

# example bin for species 2:
lines(c(rep(lenBinsStart[egBin], 2), xMasses[2]) + offset,
       c(yLengths[2], rep(lengthToMass(lenBinsStart[egBin], LWa[2], LWb[2]), 2)),
       lty=3, lwd=1, col=col2[1])
lines(c(rep(lenBinsStart[egBin+1], 2), xMasses[2]) + offset,
       c(yLengths[2], rep(lengthToMass(lenBinsStart[egBin+1], LWa[2], LWb[2]),
         2)), lty=3, lwd=1, col=col2[1])
lines(c(rep(midOfBin, 2), xMasses[2]) + offset,
       c(yLengths[2], rep(massMidOfBin[2,], 2)),
       lty=1, lwd=1, col=col2[1])

num2 = 11       # number of log2 bin breaks
log2bins = 2^(0:(num2-1))    # on unlogged axes
log2binsStart = log2bins[-num2]
log2binsEnd = log2bins[-1]   # see segments in next figure
log2binsMid = colMeans(rbind(log2binsStart, log2binsEnd))  # midpoints (unlogged)
# If want these here then copy segments code from below; probably too clutterred
#  though, and this is a slightly different issue as depends on the binning for
#  the method.
# points(rep(xMasses[1]-2, length(log2bins)), log2bins, pch="-", col="darkgrey")
# lines(rep(xMasses[1]-2, 2), c(min(log2bins), max(log2bins)), col="darkgrey")

dev.off()

# How binned body masses get assigned to logarithmic bins.
postscript("binLWlog.eps", height = 6, width = 7.5,
            horizontal=FALSE,  paper="special")

par(xaxs="i", yaxs="i")
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
lThick = 3                 # Thickness for curves
par(lend="butt")           # To have butted line caps, need for thick lines.

# Axes ranges not automated - should do if want to start changing them.
plot(100,100, xlab="", ylab="Body mass, g", xlim=c(0, 12),
     ylim=c(-60, 800), xaxt="n")      # dummy points to set up axes.

#legend("topleft", legend=c(sp2, sp1), lty=1, lwd=lThick, col=c(col2[1], col1[1]),
#       bty="n", inset=c(0.1,-0.02))
axis(2, at = seq(0, 800, by=100), labels=rep("", 9), tck=-0.02)
axis(2, at = seq(0, 800, by=50), labels=rep("", 17), tck=-0.01)
axis(1, at = seq(1, 11, by=2), labels=c("Bin 1", "Bin 2", "Bin 3", "Bin 4",
                                   "Bin 5", "Bin 6"), tck=0, col.axis="blue")

# Should really have done throughout, as clearer. Do for one species here
#  then expand to do for two species (by having columns spName, LWa, LWb).
binsSp2 = data.frame(lenStart = lenBinsStart, lenEnd = lenBinsEnd)
binsSp2 = mutate(binsSp2, lenMid = (lenStart +lenEnd)/2)
binsSp2 = mutate(binsSp2, massStart=lengthToMass(lenStart, LWa[2], LWb[2]),
     massEnd=lengthToMass(lenEnd, LWa[2], LWb[2]),
     massMid = (massStart + massEnd) / 2,
     massOfLenMid = lengthToMass(lenMid, LWa[2], LWb[2]))
        # massOfLenMid is the mass corresponding to the converted midpoint
        #  of length bins, which is actually what gets used

# Now assign each original length bin a log2 mass bin:
ind=vector()             # index of which log2 bin each mass bin ends up in
for(ii in 1:dim(binsSp2)[1])
    {
        ind[ii] = which(binsSp2[ii,"massOfLenMid"] >= log2binsStart &
               binsSp2[ii,"massOfLenMid"] < log2binsEnd)
    }
binsSp2 = cbind(binsSp2, log2binStart = log2binsStart[ind],
    log2binEnd = log2binsEnd[ind], log2binMid = log2binsMid[ind])

# Show where example bin ends up, prob just do for species 2 since clearer:

# 6 bins, therefore contain each pair in a span of 2 wide.
xVals = 1.5 + seq(0, 10, 2)     # xvals for vertical mass bins
for(egBin2 in 1:6)            # example bin number to highlight here
    {
    xVal = xVals[egBin2]                 # where to have vertical bars
    turn= xVal - 0.5       # where to turn the line - prob wrong
    # turn = 0.4*(xVals-xLog2) + xLog2
    end = xVal - 0.8                           # where to end the line
    # log2 bins:
    xLog2 = end-0.2           # where to place log2 bins
    abline(v=seq(0, 12, 2), col="grey")
    segments(x0=rep(xLog2, num2), y0=log2binsStart, y1=log2binsEnd,
         col=c("black", "grey"), lwd=thick )
    points(rep(xLog2, num2-1), log2binsMid, pch="-", cex=2, col="purple")

    segments(x0=xVal, y0=massBinsStart[2,],
         y1=massBinsEnd[2,], col=colBins[2,], lwd=thick)
    lines(c(xVal, turn, end), c(rep(binsSp2[egBin2, "massOfLenMid"],2),
       binsSp2[egBin2, "log2binMid"]),  lty=1, lwd=1, col=col2[1])
    lines(c(xVal, turn, end), c(rep(binsSp2[egBin2, "massStart"], 2),
       binsSp2[egBin2, "log2binMid"]),  lty=2, lwd=1, col=col2[1])
    lines(c(xVal, turn, end), c(rep(binsSp2[egBin2, "massEnd"], 2),
       binsSp2[egBin2, "log2binMid"]),  lty=2, lwd=1, col=col2[1])
    # Arrow for the final part:
    Arrows(x0=turn, y0=binsSp2[egBin2, "massOfLenMid"], x1 = end,
       y1 = binsSp2[egBin2, "log2binMid"], arr.adj=1,
       col = col2[1], lwd=1, arr.lwd=0.1, arr.type="triangle")
   }
dev.off()

# save.image(file="binLW.RData")   # save everything
