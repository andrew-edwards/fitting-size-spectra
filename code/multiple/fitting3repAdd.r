# fitting3repAdd.r - results from fitting3rep.r, and adding gold histograms
#  for results from a second set of simulations, as documented
#  in READdataCode.md.
#  Need to re-run fitting3rep.r if want new results. 23rd November 2015.
# fitting3rep.r - as for fitting2rep.r, but now plotting histograms in terms
#  of estimates of b rather than estimates of slope, since they will be clearer.
#  13th July 2015.



# Saving resulting figure as fitting3repAdd**.eps, where ** is:

#     **          xmax1    seed1   xmax2   seed2
#     a            10^3     42      10^4    42
#     b            10^3     42      10^4    43
#     c            10^3     42      10^3    43
#     d            10^4     42      10^4    43
# First results to load in are xmax2, seed2 (resaved as .set2),
#  because they will
#  be plotted first (in the background), and wanted xmax1, seed1 to
#  be prominent (and so the rest of the .RData variables refer to that
#  main one).
# So change the two load() commands and the .eps filename.

rm(list=ls())
# xmax2 and seed2:
load("xmax10000/fitting3rep10000.RData")    # a
# load("xmax10000b/fitting3rep10000b.RData")    # b, d   (10000b is seed=43)
# load("xmax1000-43/fitting3rep1000-43.RData")    # c

# Save results as ...set2
toRemove = ls()
Llin.rep.set2 = Llin.rep
LT.rep.set2 = LT.rep
LTplus1.rep.set2 = LTplus1.rep
LBmiz.rep.set2 = LBmiz.rep
LBbiom.rep.set2 = LBbiom.rep
LBNbiom.rep.set2 = LBNbiom.rep
LCD.rep.set2 = LCD.rep
MLE.rep.set2 = MLE.rep
rm(list = toRemove)    # Remove everything loaded in, since
                       #  variable names are repeated in fitting3rep.Rdata.

# xmax1 and seed1:
load("fitting3rep.RData")     # a, b, c
# load("xmax10000/fitting3rep10000.RData")    # d
source("../PLBfunctions.r")

xrange = c(-3.5, 0.5)         # range of x-axis for histograms - actually to define bins
# xbigticks = -c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25)    # Change to depend on brange
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
postscript("fitting3repAdda.eps", height = figheight, width = figwidth,
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
xLim = c(-3.2, -1.5)
xLimLlin = xLim + 2   
ylimA = c(0, 5500)
ylimLlin = c(0, 10000)
# ylimB = c(0, 3250)
cexAxis = 0.9        # font size for axes labels to make the y ones fit okay

bordCol1 = "blue"  # colour for border (no shading) for set1
shadCol= "gold"   # shading colour for the set2 histograms. Gold is good.
                 #  rgb(1, 0, 0, 0.5) for transparent blue doesn't work on .eps
bordCol2 = shadCol    # border colour for the set2 histograms

# Llin.rep has different breakhist so that 0 is a breakpoint. Same width as
#  others, just shifted.
# Need to do set2 histograms  first, to then get overdrawn
#  by the more prominent set1 (usually main xmax=1000) results.
hist(Llin.rep.set2, xlim=xLimLlin, breaks=breakshist - breakshist[31],
  xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimLlin, col=shadCol, border=bordCol2)  #  ylim=c(0,1100),
hist(Llin.rep, breaks=breakshist - breakshist[31], add=TRUE, border=bordCol1)
axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks
axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)
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


hist(LT.rep.set2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)  #  ylim=c(0,1100), 
hist(LT.rep, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(b) LT", bty="n", inset=inset)

hist(LTplus1.rep.set2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)  #  ylim=c(0,1100), 
hist(LTplus1.rep, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(c) LTplus1", bty="n", inset=inset)

hist(LBmiz.rep.set2 -1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)  #  ylim=c(0,1100),
hist(LBmiz.rep -1, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(d) LBmiz", bty="n", inset=inset)

hist(LBbiom.rep.set2 -2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)  #  ylim=c(0,1100),
hist(LBbiom.rep -2, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(e) LBbiom", bty="n", inset=inset)

hist(LBNbiom.rep.set2 -1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)  #  ylim=c(0,1100),
hist(LBNbiom.rep -1, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(f) LBNbiom", bty="n", inset=inset)


hist(LCD.rep.set2 -1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)
hist(LCD.rep -1, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(g) LCD", bty="n", inset=inset)

hist(MLE.rep.set2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA, col=shadCol, border=bordCol2)
hist(MLE.rep, breaks=breakshist, add=TRUE, border=bordCol1)
histAxes()
# inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(h) MLE", bty="n", inset=inset)

# text(-2.5, 0, expression(paste("Estimate of ", italic(b))))
mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)
# mtext("hello", side=1, outer=TRUE, line=-1)

dev.off()

# save.image(file = "fitting3repAdd.RData")
