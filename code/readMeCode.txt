readMeCode.txt - readMe file for R code for manuscript
 'Testing and recommending methods for fitting size spectra to data'
 by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum
 and Julia L. Blanchard, submitted to Method in Ecology and Evolution. 


**The code that fits all the methods (e.g. fitting3rep.r) can take a short
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

fitting3rep.RData - results from fitting3rep.r, to save having to re-run it.

fitting3repAdd.r - constructing Figure 3, combining simulation results from
 two sets of 10,000 simulated data sets.

fitting3conf.r - Figure 4 plots of confidence intervals, and Figure A.4 for
 MLEfix method.

fitting3bmaxx.r - Figure A.2, showing relationship between MLE of b and MLE 
 of x_max for the 10,000 simulated data sets from Figure 3. And Figure A.5 for
 MLEfix method.

The following are the sensitivity analyses. Some are just modifications 
 of the above code, with new a value of a parameter (e.g. xmax)
 and renaming of anything that is saved (.RData and .eps files), and 
 changing of axes sizes where necessary (and maybe some other necessary tweaks).
 The above code may have been updated (written more clearly) after some of
 the following were done.

code/multiple/xmax10000/  - simulate 10,000 datasets with xmax = 10,000
************************

fitting2-10000new.r - as for fitting2.r but for xmax=10,000, to produce 
 Figures A.6 and A.7.

fitting3rep10000.r - as for fitting3rep.r but for xmax=10,000, to give
 gold histograms in Figure 3 and results in Table A.1.

fitting3rep10000.RData - results from fitting3rep10000.r to save having
 to re-run it.

fitting3conf10000.r - as for fitting3conf.r but for xmax=10,000, to give
 Figure A.8.

fitting3bmaxx10000.r - Figure A.9 for the MLE and MLEfix methods, with 
 xmax = 10,000.


code/bMinus25/    b = -2.5
**************

fitting2bMinus25.r - Figure A.10.

fitting3rep-bMinus25.r - Figure A.11 and Table A.2.

fitting3rep-bMinus25.RData - results from fitting3rep-bMinus25.r.

fitting3conf-bMinus25.r - Figure A.12

code/bMinus15/   b = -1.5
**************

fitting2bMinus15.r - Figure A.13.

fitting3rep-bMinus15.r - Figure A.14 and Table A.3.

fitting3rep-bMinus15.RData - results from fitting3rep-bMinus15.r.

fitting3conf-bMinus15.r - Figure A.15


code/bMinus05/   b = -0.5
**************

fitting2bMinus05.r - Figure A.16.

fitting3rep-bMinus05.r - Figure A.17 and Table A.4.

fitting3rep-bMinus05.RData - results from fitting3rep-bMinus05.r.

fitting3conf-bMinus05.r - Figure A.18

code/n10000/     n = 10000
************

Redoing main results with sample size (number of individual measurements) 
 increased to 10000.

fitting2n10000.r - Figure A.19

fitting3rep-n10000.r - Figure A.20 and Table A.5

fitting3rep-n10000.RData - results from fitting3rep-n10000.r.

fitting3conf-n10000.r - Figure A.21


code/MLEbin/    MLEbin method for likelihood when the data are already binned
************

fitting3repMLEbin.r - same simulated data sets as in fitting3rep.r, but
 binning the data and then applying likelihood.

fitting3repMLEbin.Rdata - results from fitting3repMLEbin.f

fitting3confMLEbin.r  - confidence intervals for MLEbin method, to give Figure 5.


code/recommend/
***************

recommend.r - Figure 6, recommended MLE calculations and resulting plots of 
 data and fitted size spectrum.

--

To test sensitivity of results to xmax, b etc. (as I've done above), create a
 new folder, and then copy in, rename and edit the files fitting2.r, 
 fitting3rep.r and fitting3conf.r (or ones above that are closer to what
 you are testing - e.g. use b=-2.5 if testing b=-2.6). 
Edits include: change required parameter value, change .eps and .RData
 filenames, and then may need to manually edit axes since it's hard to
 fully automate them, particularly the barplot with a gap, as in the 
 function  gap.barplot.cust.

To re-run code with a different seed, say 43, just change set.seed(42)
 to set.seed(43) and make sure redo.simulation=TRUE (if it's there) so
 that it doesn't just load in an already saved .RData file. It's best to first
 move the code to a new directory.
 


















 

