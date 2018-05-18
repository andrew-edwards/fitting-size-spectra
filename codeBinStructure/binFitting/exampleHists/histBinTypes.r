# histBinTypes.r - histograms of a single data set using the different
#  binTypes. 15th June 2016.

# fitting2.r - replacing the fitting methods with functions, that are all
#  in PLBfunctions.r. Some axes limits may be chosen explicitly for the
#  simulated data set and so are not automated. 5th June 2015.

rm(list=ls())
print(date())

# require(Hmisc)            # for subplot()
source("../../../code/PLBfunctions.r")
source("../../countsFunctions.r")

n = 1000                  # sample size
b.known = -2              # known fixed value of b
  xmin.known = 1            # known fixed value of xmin. Make it a power of 2
  if(!is.wholenumber(log2(xmin.known)) | !is.wholenumber(xmin.known))
       { stop("If want xmin.known to not be an integer and not be an integer
                 power of 2 then need to edit binData; may just need to think
                 about whether can just set startInteger = FALSE, but will
                 need to edit binData to be able to define the binWidth for
                 2k method. Will need to think about this for real data; maybe
                 best to just remove the first bin (and a fit a range that is
                 encompassed by the binned data).")
      }
xmax.known = 1000         # known fixed value of xmax

set.seed(42)              # Same seed as for original simulations in
                            #  first manuscript.

binType = list(1, 5, 10, "2k", "2k", 50)
                                   # Must be numeric, for linear bins, or "2k"
                                   #  so define as a list
freq = rep(TRUE, length(binType))  # Whether to plot frequency or density
freq[max(which(binType == "2k"))] = FALSE # Plot density for second 2k

yMax = rep(xmax.known, length(binType))  # Same y-axes
yMax[binType=="1"] = 600                 # Shorter y-axis
yMax[binType=="2k"] = 600
yMax[!freq] = 0.6                        # For 2k when plotting density

binType.name = binType
binType.name[!binType %in% "2k"] = paste0("Linear ",
                  binType[which(!binType %in% "2k")])

x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)

figheight = 7
figwidth = 5.7    # 5.7 inches for JAE

par(mfrow=c(1,2))
for(j in 1:length(binType))                    # Loop over binning type
  {
    bins.list = binData(x, binWidth=binType[[j]])  # 1, 2, 5 or "2k", etc.
    num.bins = dim(bins.list$binVals)[1]   # may not need

    binBreaks = bins.list$binVals[,"binMin"]$binMin   # Loses the column names
    maxOfMaxBin = bins.list$binVals[num.bins, "binMax"]$binMax
    binBreaks = c(binBreaks, maxOfMaxBin) # Append endpoint of final bin

    binCounts = bins.list$binVals[,"binCount"]$binCount
    # binMids = bins.list$binVals[,"binMid"]$binMid  # Midpoints of bins
    figname = paste0("histBinType", binType[j])
    if(!freq[j]) figname = paste0(figname, "density")
    figname = paste0(figname, ".eps")

    postscript(figname,
           height = figheight, width = figwidth,
           horizontal=FALSE, paper="special")
    h = hist(x, breaks=binBreaks, # xaxs="i", yaxs="i",
        xlab=expression(paste("Values, ", italic(x))),
        ylim=c(0,yMax[[j]]),
        border=NA,        # set border or lty
        # lty="blank",
        col="red", lwd=1,
        main=binType.name[j],
        freq=freq[j])              # Will give warning for 2k
        # col=rep("red", length(binBreaks)-1)
    if(binType[j] == "FALSE") # "2k")
       {
       subplot( hist(x, breaks=binBreaks,
            border=NA,                # or lty="blank",
            col="red", lwd=1, main = "", xlab="", ylab="",
            ylim=c(0,0.125*yMax[j]),
            freq=freq[j]),
            350, 0.75*yMax[j], size=c(2,2))
       print(binCounts)
       print(binBreaks)
       print(diff(binBreaks))
       print(binCounts/diff(binBreaks))
       aa = binCounts/diff(binBreaks)
       bbbb = h
       }
    dev.off()
    if(max(abs(h$counts - binCounts)) > 0) stop("Check h$counts.")
  }   # End of for(j in 1:binTypes) loop


