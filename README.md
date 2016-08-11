# fitting-size-spectra

This repository contains R code for the paper 

**Testing and recommending methods for fitting size spectra to data** 

by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard

that has been accepted for publication in *Methods in Ecology and Evolution*. 

Please email me if you would like a copy of the paper, and please cite it if you end up using the code.

The aim of sharing this code is so that others can repeat (and extend) our simulation study, and also analyse their own data.

To download the code from the GitHub site just click the 'Clone or Download' button (near the top on the right) and select 'Download ZIP'. If you use GitHub then feel free to fork and even adapt the code.

To *exactly* reproduce the results in the paper
download release version 1.0.0, as finalised on 24th May 2016, and submitted to *Methods in Ecology and Evolution* with the revised version of the manuscript. [Click on 'release' tab in the GitHub site]. Later updates of the code will have some generalisation in functions that should not affect the older code (but I just won't re-test it all.

There are also functions (in **code/PLBfunctions.r**) that may be of more general use, such as **logTicks()** for adding tick marks to a log-log plot, and **legJust()** for right-justifying a legend (functionalising an example in ?legend). 

Andrew Edwards. 

Andrew.Edwards@dfo-mpo.gc.ca (may change later in 2016).

http://www.chebucto.ns.ca/~english 

# Repository Contents

**README.md** - this file.

**.gitignore** - ignore this if you don't know what it is.

**code/** - directory containing all the R code.

**code/readMeCode.txt** - readme file for the code directory. Contains instructions and details of the files and further subdirectories.

The subdirectories of **code/** are summarised below, but see **readMeCode.txt** for full details.

**code/single/** - testing the eight methods on a single data set.

**code/multiple/** - testing the eight methods on multiple (10,000) simulated data sets. Contains subdirectories for the sensitivity tests.

**code/MLEbin/** - MLEbin method for likelihood when the data are only available in binned form.

**code/recommend/** - recommended likelihood calculations and resulting plots of data and fitted size spectrum (Figure 6).

# Exactly reproducing results

To ***exactly*** reproduce the results in our manuscript use release version 1.0.0 (see above). You need to run each R code (that generates random numbers) in a fresh R window. But why, when the seed is set in each .r file?

**TO FINISH**