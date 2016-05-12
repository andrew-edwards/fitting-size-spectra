# fitting3bmaxx10000.r - as for fitting3bmaxx.r but with xmax=10000. 17/7/15.

# fitting3bmaxx.r - plotting MLE for b versus xmax (=max(x) for each (or
#  maybe a subset) of the 10,000 simulated data sets. For original MLE
#  and for MLEfix, which uses xmax=xmax.known for each dataset. 17th July 2015

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

load("fitting3rep10000.RData")
source("../../PLBfunctions.r")  # to load in required functions (probability
                          #  functions for PL and PLB and more, including
                          #  gap.barplot.andy)

# To give same ranges on y axes: 
# range(MLE.rep)      -2.123247 -1.879260
# range(MLEfix.rep)   -2.127078 -1.885421
yLim = c(-2.13, -1.87)

postscript("fitting3bmaxx10000.eps",  height = 5.4, width = 5.7,
           horizontal=FALSE,  paper="special")  

# par(mai=c(0.6, 0.6, 0.15, 0.3))
# mgpVals = c(2, 0.5, 0)            # mgp values   2.0, 0.5, 0

lm.MLE = lm(MLE.rep ~ MLE.rep.xmax)

plot(MLE.rep.xmax, MLE.rep, pch=20, cex=0.3, xlab="MLE of xmax",
     ylab="MLE of b", xlim=c(0, xmax.known), ylim=yLim)
lm.line(x.vector=range(MLE.rep.xmax), lm.results=lm.MLE, col="red")

# Based on RBR, plot the confidence intervals, but they're quite narrow:
xmax.inc = seq(min(MLE.rep.xmax), max(MLE.rep.xmax), length=100)
p.conf = predict(lm.MLE, newdata=data.frame(MLE.rep.xmax=xmax.inc),
    interval="confidence")
matlines(xmax.inc, p.conf[ ,c("lwr", "upr")], col="red", lty=2)
## pVal = summary(lm.MLE)$coeff["MLE.rep.xmax",4]
## if (pVal[i] > 0.05) regCol = "grey" else regCol="red"
inset = c(-0.05, -0.02) #c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(a) MLE method", bty="n", inset=inset)

dev.off()

print(paste("Correlation coefficient is", cor(MLE.rep.xmax, MLE.rep)))
print(summary(lm.MLE))      

postscript("fitting3bmaxxMLEfix10000.eps",  height = 5.4, width = 5.7,
           horizontal=FALSE,  paper="special")  
# par(mai=c(0.6, 0.6, 0.15, 0.3))
# mgpVals = c(2, 0.5, 0)            # mgp values   2.0, 0.5, 0

lm.MLEfix = lm(MLEfix.rep ~ MLE.rep.xmax)   #MLE.rep.xmax is max(x) for each

plot(MLE.rep.xmax, MLEfix.rep, pch=20, cex=0.3, xlab="Max. of data set",
     ylab="MLE of b", xlim=c(0, xmax.known), ylim=yLim)
lm.line(x.vector=range(MLE.rep.xmax), lm.results=lm.MLEfix, col="red")

# Based on RBR, plot the confidence intervals, but they're quite narrow:
# xmax.inc = seq(min(MLE.rep.xmax), max(MLE.rep.xmax), length=100)
p.confFix = predict(lm.MLEfix, newdata=data.frame(MLE.rep.xmax=xmax.inc),
    interval="confidence")
matlines(xmax.inc, p.confFix[ ,c("lwr", "upr")], col="red", lty=2)
## pVal = summary(lm.MLEfix)$coeff["MLE.rep.xmax",4]
## if (pVal[i] > 0.05) regCol = "grey" else regCol="red"
legend("topleft", "(b) MLEfix method", bty="n", inset=inset)
dev.off()

print(paste("For MLEfix, correlation coefficient is", cor(MLE.rep.xmax, MLEfix.rep)))
print(summary(lm.MLEfix))      


