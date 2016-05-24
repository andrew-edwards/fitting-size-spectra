# recommend.r - recommended way of estimating exponent b and of plotting
#  the results. Uses the same simulated data set as start of paper
#  as an example. 23rd November 2015.

rm(list=ls())
print(date())

n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 1000         # known fixed value of xmax

source("../PLBfunctions.r")
                          # to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy)

# Sample from PLB distribution:
set.seed(42)      # To get the same observations for each run of code.
                  # 8 bins, only up to 400 with 2 empty.
newdata = TRUE    # TRUE - generate new data, FALSE - load in previous data
                  #  (especially if the calculations take a while but you only
                  #  want to tweak the figures).
                  #  Do NOT change seed and set to TRUE without editing
                  #  filenames for saving results and figures (below).
if(newdata)
  {
  x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)
  } else
  {
  load(file="recommend.RData")    # or load in data for x. Replace this with
                                  #  your own data set.
  }


# x is a vector of individual fish sizes (here it is body masses because
#  we then calculate the biomass size spectrum)

log.x = log(x)                      # to avoid keep calculating
sum.log.x = sum( log.x ) 
xmin = min(x)
xmax = max(x)

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

bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.00001) 
    
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

figheight = 4.2 # 7 for 4x2 figure
figwidth = 2.85

postscript("recommend.eps", height = figheight, width = figwidth,
           horizontal=FALSE, paper="special")
par(mfrow=c(2,1))
oldmai = par("mai")    #  0.95625 0.76875 0.76875 0.39375  inches I think,
                       #   think may be indpt of fig size
par(mai=c(0.4, 0.5, 0.05, 0.3), cex=0.7)  # Affects all figures if don't change again
mgpVals = c(1.6,0.5,0)            # mgp values   2.0, 0.5, 0

# Notation:
# hAAA - h(istrogram) for method AAA.

inset = c(0, -0.04)     # inset distance of legend

# LBNbiom method - on biomass, not counts. Calculates
#  log2 bins of bodymass, sum the total biomass in each bin, normalises
#  biomasses by binwidths, fits regression to log10(normalised biomass) v
#  log10(midpoint of bin), but we're not going to use the regression slope
#  here, but we need the binning.
hLBNbiom.list = LBNbiom.method(x)

# Plotting Figure 6(a)

plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals, xlim=c(0,2.7),
     ylim=c(0,3), yaxt="n")

if(min(hLBNbiom.list[["binVals"]]$log10binMid) < par("usr")[1]
    | max(hLBNbiom.list[["binVals"]]$log10binMid) > par("usr")[2])
   { stop("fix xlim for LBNbiom method")} 
if(min(hLBNbiom.list[["binVals"]]$log10totalBiomNorm) < par("usr")[3]
    | max(hLBNbiom.list[["binVals"]]$log10totalBiomNorm) > par("usr")[4])
   { stop("fix ylim for LBNbiom method")} 

axis(2, at = 0:3, mgp=mgpVals)
axis(2, at = c(0.5, 1.5, 2.5), mgp=mgpVals, tcl=-0.2, labels=rep("", 3))

x.PLB = seq(min(x), max(x), length=1000)     # x values to plot PLB. Note
                           # that these encompass the data, and are not based
                           # on the binning (in Figure 6 the line starts as
                           # min(x), not the first bin.

B.PLB = dPLB(x.PLB, b = PLB.bMLE, xmin=min(x.PLB),
    xmax=max(x.PLB)) * length(x) * x.PLB
       # The biomass density, from equation (7), using the MLE for b.

lines(log10(x.PLB), log10(B.PLB), col="red")

# To see all curves for b within the confidence interval:
#for(i in 1:length(bIn95))  
#    {
#        lines(log10(x.PLB), log10( dPLB(x.PLB, b = bIn95[i], xmin=min(x.PLB),
#            xmax=max(x.PLB)) * length(x) * x.PLB), col="blue", lty=2)
#    }

# To add just the curves at the limits of the 95% confidence interval of b:
for(i in c(1, length(bIn95))) 
    {
        lines(log10(x.PLB), log10( dPLB(x.PLB, b = bIn95[i], xmin=min(x.PLB),
            xmax=max(x.PLB)) * length(x) * x.PLB), col="red", lty=2)
    }

# This would plot the estimate from the LBNbiom method, as in Figure 2(f).
# lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])

legend("topright", "(a)", bty="n", inset=inset)

# Plotting Figure 6(b):
plot(sort(x, decreasing=TRUE), 1:length(x), log="xy",
     xlab=expression(paste("Values, ", italic(x))),
     ylab=expression( paste("Number of ", values >= x)), mgp=mgpVals,
     xlim = c(xmin, xmax), ylim = c(1, n), axes=FALSE)
xLim = 10^par("usr")[1:2]
yLim = 10^par("usr")[3:4]

logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500))   # Tick marks.

y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                  xmax = max(x.PLB))) * length(x)    
lines(x.PLB, y.PLB, col="red") #, lty=5)

# legend("topright", paste("(b) Exponent=", signif(PLB.bMLE, 3)), bty="n",
#       inset=inset)
legend("topright", "(b)", bty="n", inset=inset)


# To see all curves for b within the confidence interval:
#for(i in 1:length(bIn95))  
#    {
#      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
#                  xmax = max(x.PLB))) * length(x), col="blue", lty=2)
#    }  

# To add just the curves at the limits of the 95% confidence interval of b:
for(i in c(1, length(bIn95)))
    {
      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
                  xmax = max(x.PLB))) * length(x), col="red", lty=2)
    }  

dev.off()

print(paste("MLE of b is", round(PLB.bMLE, dig=5)))
print(paste("95% confidence interval of b is (", round(PLB.MLE.bConf[1], dig=5),
            ", ", round(PLB.MLE.bConf[2], dig=5), ")", sep=""))

save.image(file = "recommend.RData")


