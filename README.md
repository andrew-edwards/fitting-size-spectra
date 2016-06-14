# fitting-size-spectra

***generalising*** branch - updating some code for second manuscript, that should still work fine with original code on ***master***. If you are downloading code, please go to the ***master*** branch.



This repository contains R code for the manuscript 

**Testing and recommending methods for fitting size spectra to data** 

by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard

that has been submitted to *Methods in Ecology and Evolution*. 

Please email me if you would like a copy of the manuscript, and please cite the manuscript if you end up using the code.

The aim of sharing this code is so that others can repeat (and extend) our simulation study, and also analyse their own data.

To download the code from the GitHub site just click the 'Download ZIP' button (near the top on the right). If you use GitHub then feel free to fork and even adapt the code.

This is release version 1.0.0, as finalised on 24th May 2016, and submitted to *Methods in Ecology and Evolution* with the revised version of the manuscript.

Andrew Edwards. 

Andrew.Edwards@dfo-mpo.gc.ca (will change later in 2016).

www.chebucto.ns.ca/~english 

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
