readMeCode.txt - readMe file for R code for manuscript
 'Testing and recommending methods for fitting size spectra to data'
 by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum
 and Julia L. Blanchard, submitted to Method in Ecology and Evolution. 


At some point, edit:

Note to users:

The code that fits all the methods (e.g. fitting3rep.r) can take a short
 while to run, and so I have included the .RData files for those.

I have tried to keep required packages to a minimum, but you will need:
plotrix
dplyr
(and their dependencies).



code/
*****

readMeCode.txt - this file

PLBfunctions.r - collection of documented functions called by the other 
 R code.

code/single/ - simulate a single data set and fit spectra using the eight methods
************

fitting2.r - simulates a data set and then fits spectra using eight methods,
 producing Figures 1, 2 and A.1.



code/multiple/ - simulate 10,000 data sets and fit using the eight methods
**************

fitting3rep.r - results from 10,000 simulated data sets, to give the blue
 histograms in Figure 3 and the main results in Table 2. Also does the
 MLEfix method and plots Figure A.2 (or A.3?? - **Check).




fitting3conf.r - Figure 4 plots of confidence intervals, and Figure A.4 for
 MLEfix method.








Below here is from readMEErevised.txt in size-spectra-methods. When I move something to fitting-size-spectra I will then move it to above here and document it in the directory structure.

Contents of MEEsubmit.zip
*************************



fitting3rep.RData - results from fitting3rep.r, to save having to re-run it.

fitting3repAdda.r - constructing Figure 3, combining simulation results from
 two sets of 10,000 simulated data sets.


recommend.r - Figure 6, recommended presentation of data and fitted 
 size spectrum.

fitting3bmaxx.r - Figure A.1, showing relationship between MLE of b and MLE 
 of x_max for the 10,000 simulated data sets from Figure 3. And Figure A.4 for
 MLEfix method.

xmax10000\  directory: xmax = 10,000
************************************

Some of these are just modifications of the above code, with new value of xmax
 and renaming of anything that is saved (.RData and .eps files), and changing of
 axes sizes where necessary (and maybe some other necessary tweaks).


fitting2-10000new.r - as for fitting2.r but for xmax=10,000, to produce 
 Figures A.5 and A.6.

fitting3rep10000.r - as for fitting3rep.r but for xmax=10,000, to give
 gold histograms in Figure 3 and results in Table A.1.

fitting3rep10000.RData - results from fitting3rep10000.r to save having
 to re-run it.

fitting3conf10000.r - as for fitting3conf.r but for xmax=10,000, to give
 Figure A.7.

fitting3bmaxx10000.r - Figure A.8 for the MLE and MLEfix methods, with 
 xmax = 10,000.


bMinus25\  directory: b = -2.5
******************************

fitting2bMinus25.r - Figure A.9.

fitting3rep-bMinus25.r - Figure A.10 and Table A.2.

fitting3rep-bMinus25.RData - results from fitting3rep-bMinus25.r.

fitting3conf-bMinus25.r - Figure A.11


bMinus15\  directory: b = -1.5
******************************

fitting2bMinus15.r - Figure A.12.

fitting3rep-bMinus15.r - Figure A.13 and Table A.3.

fitting3rep-bMinus15.RData - results from fitting3rep-bMinus15.r.

fitting3conf-bMinus15.r - Figure A.14


**TO ADD TO .7z, below here*****:

bMinus05\  directory: b = -0.5
******************************

fitting2bMinus05.r - Figure A.**.

fitting3rep-bMinus05.r - Figure A.** and Table A.**.

fitting3rep-bMinus05.RData - results from fitting3rep-bMinus05.r.

fitting3conf-bMinus05.r - Figure A.**



n10000\  directory: n = 10000
*****************************

Redoing main results with sample size (number of individual measurements) 
 increased to 10000.

fitting2n10000.r - Figure A.**

fitting3rep-n10000.r - Figure A.** and Table A.**

fitting3rep-n10000.RData - results from fitting3rep-n10000.r.

fitting3conf-n10000.r - Figure A.**


fitting3repMLEbin.r - for Figure **
fitting3repMLEbin.Rdata 
fitting3confMLEbin.r 


To re-run code with a different seed, say 43, just change set.seed(42)
 to set.seed(43) and make sure redo.simulation=TRUE (if it's there) so
 that it doesn't just load in an already saved .RData file.







 

