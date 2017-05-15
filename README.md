# fitting-size-spectra

This repository contains R code for the paper 

**Testing and recommending methods for fitting size spectra to data** 

by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard

published in ***Methods in Ecology and Evolution*** (2017, 8:57-67). The paper is freely available at

<http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full>

The aim of sharing this code is so that others can repeat (and extend) our simulation study, and also analyse their own data. I have endeavoured to properly document the code.

To download the code from the GitHub site just click the 'Clone or Download' button (near the top on the right) and select 'Download ZIP'. Make a note of the 'Latest commit' number in case you have any questions for me. If you use GitHub then feel free to fork and adapt the code.

If you want to run new simulations or apply the code to your own data then just use the latest version that will be automatically displayed on the GitHub site.

To *exactly* reproduce the results in the paper
you should download release version 1.0.0 (click on 'release' tab in the GitHub site) and see the notes in **code/readReRun.txt** concerning the seed and how **require(dplyr)** generates a random number. Later updates of the code will have some generalisation of functions that should not affect the older code (but I won't re-test it all for back compatibility). 

There are functions (in **code/PLBfunctions.r**) that may be of more general use, such as **logTicks()** for adding tick marks to a log-log plot, and **legJust()** for right-justifying a legend (based on an example in **?legend**). 

If you have problems with the code then please contact me. Some of it has been independently used by co-author James Robinson, who got it working fine for his own data.

Thanks,

Andrew Edwards. 

<http://www.chebucto.ns.ca/~english>

Andrew.Edwards@dfo-mpo.gc.ca - this may change at some point, and if it does not work then try andrew.edwards.dfo@gmail.com which will always automatically forward to my correct email.

# Repository Contents

**README.md** - this file (can be read in a Markdown viewer or simply any text editor).

**.gitignore** - ignore this if you don't know what it is.

**code/** - directory containing all the R code.

**code/readMeCode.txt** - readme file for the code directory. Contains instructions and details of the files and further subdirectories.

The subdirectories of **code/** are summarised below, but see **readMeCode.txt** for full details.

**code/single/** - testing the eight methods on a single data set.

**code/multiple/** - testing the eight methods on multiple (10,000) simulated data sets. Contains subdirectories for the sensitivity tests.

**code/MLEbin/** - MLEbin method for likelihood when the data are only available in binned form.

**code/recommend/** - recommended likelihood calculations and resulting plots of data and fitted size spectrum (Figure 6).

15th May 2017 - corrected a minor error regarding the likelihood function for the MLEbin method when b=-1. In practice this should not have any effect. See Issue #7. 