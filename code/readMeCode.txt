readMeCode.txt - readMe file for R code for manuscript
 'Testing and recommending methods for fitting size spectra to data'
 by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum
 and Julia L. Blanchard, submitted to Method in Ecology and Evolution. 

 All code by Andrew Edwards, and available at 
  https://github.com/andrew-edwards/fitting-size-spectra .
  andrew.edwards@dfo-mpo.gc.ca (though email will change soon so please check
  website at  www.chebucto.ns.ca/~english ).

This is version 1.0.0, as submitted to the journal with the revised version
 of our manuscript. I will try and label this version as a release on
 the GitHub site. 24th May 2016.

I have tried to keep required packages to a minimum, but you will need:
  plotrix  (only for histograms with gaps, as for Figure 1)
  dplyr
  (and their dependencies).

The resulting figures (as .eps postscript files) are also included so that
 code can be independnetly re-run and the results easily compared with my 
 original figures. I have included the .RData file for the results for 
 the main simulation of 10,000 data sets (multiple/fitting3rep.r) 
 and for the xmax=10000 simulation (multiple/xmax10000/fitting3rep10000.r)
 because the code can take a while to run and the two resulting .RData files
 are used to produce Figure 3. The remainaing .RData files are have
 not been included because they are generally large.

The main figures in the manuscript can be found in the following directories:

Figures 1 and 2 - code/single/
Figures 3 and 4 - code/multiple/
Figure 5 - code/MLEbin/
Figure 6 - code/recommend/

So to use the MLE method to analyse your own data and plot results as per
 our Figure 6, see code/recommend/ . If your data are binned then you will
 need some of code/MLEbin/ .

Code was mainly developed under R version 3.1.0, although I then upgraded to
 version 3.2.3 in January 2016, and so I  have re-run the code in 3.2.3
 to verify that  I get the same results.

I have functionalised code where practical, though I did not go back everywhere
 and replace original non-function code with functions; further improvements 
 are no doubt possible. 

The full list of directories and code is now given.

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
 MLEfix method and plots Figure A.3.

fitting3rep.RData - results from fitting3rep.r, to save having to re-run it.

fitting3repAdd.r - constructing Figure 3, combining simulation results from
 two sets of 10,000 simulated data sets (fitting3rep.r and 
 xmax10000/fitting3rep10000.r).

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

fitting3rep10000.RData - results from fitting3rep10000.r, which gets
 called in multiple/fitting3repAdd.r to produce Figure 3. 

fitting3conf10000.r - as for fitting3conf.r but for xmax=10,000, to give
 Figure A.8.

fitting3bmaxx10000.r - Figure A.9 for the MLE and MLEfix methods, with 
 xmax = 10,000.

compareXmax.r - comparing sample of 1,000 random PLB numbers
 for xmax=1,000 and xmax=10,000 with same seed, as documented in
 Section A.2.8 in the Appendix.


code/bMinus25/    b = -2.5
**************

fitting2bMinus25.r - Figure A.10.

fitting3rep-bMinus25.r - Figure A.11 and Table A.2.

fitting3conf-bMinus25.r - Figure A.12


code/bMinus15/   b = -1.5
**************

fitting2bMinus15.r - Figure A.13.

fitting3rep-bMinus15.r - Figure A.14 and Table A.3.

fitting3conf-bMinus15.r - Figure A.15


code/bMinus05/   b = -0.5
**************

fitting2bMinus05.r - Figure A.16.

fitting3rep-bMinus05.r - Figure A.17 and Table A.4.

fitting3conf-bMinus05.r - Figure A.18


code/n10000/     n = 10000
************

Redoing main results with sample size (number of individual measurements) 
 increased to 10000.

fitting2n10000.r - Figure A.19

fitting3rep-n10000.r - Figure A.20 and Table A.5

fitting3conf-n10000.r - Figure A.21


code/MLEbin/    MLEbin method for likelihood when the data are already binned
************

fitting3repMLEbin.r - same simulated data sets as in fitting3rep.r, but
 binning the data and then applying likelihood.

fitting3repMLEbin.Rdata - results from fitting3repMLEbin.r, since this file is 
 small.

fitting3confMLEbin.r  - confidence intervals for MLEbin method, to give Figure 5.

   /MLEbinRerun/ - rerunning MLEbin code after fixing b=-1 error (see Issue #7
   *************     on GitHub page). Not committing this folder to GitHub as
                     it would just add clutter and is not important. 

   /MLEbinRerun2/ - rerunning MLEbin code after completely fixing b=-1 error
   **************     in December 2012 (see Issue #7); again not committing.

   /MLEbin-bMinus1/ - testing with b=-1 after fixing the errors in
   ****************    Issue #7, but it appears to make no difference
                       as in the nlm search you still never have exactly
                       b == -1 (results summary is same when running from
                       commit 36a43ff, before fixing the Issue).

code/recommend/
***************

recommend.r - Figure 6, recommended MLE calculations and resulting plots of 
 data and fitted size spectrum.

recommend.eps - Figure 6.

--

To test sensitivity of results to xmax, b etc. (as I've done above), create a
 new folder, and then copy in, rename and edit the files fitting2.r, 
 fitting3rep.r and fitting3conf.r (or ones above that are closer to what
 you are testing - e.g. use b=-2.5 if testing b=-2.6). 

Such edits include: change required parameter value, change .eps and .RData
 filenames, and then may need to manually edit axes since it's hard to
 fully automate them, particularly the barplot with a gap, as in the 
 function  gap.barplot.cust.

To re-run code with a different seed, say 43, just change set.seed(42)
 to set.seed(43) and make sure redo.simulation=TRUE (if it's there) so
 that it doesn't just load in an already saved .RData file. It's best to first
 move the code to a new directory.


Most files have a  
   redo.simulation = TRUE   
or
   redo.simulation = FALSE
option at the start, and they will be mostly set to FALSE to just load in
the simulation results, because I would have been tweaking the figures for
publication. So obviously set to TRUE for the first run, until you have an
.RData file that can then be loaded in. 

I have re-run all code in R version 3.2.3 to ensure results are as stated in the
 manuscript. I strangely find that if I re-run code such as fitting3rep.r
 (that generates 10,000 sets of 1,000 random numbers) from a new R console,
 I get a different final set of random numbers answer than re-running in a
 console in which I had just run fitting2.r. The first set of 1,000 random 
 numbers is the same, but the last is shifted along by one. I will try and
 create a minimum working example to investigate this.

This occurs even though I have
  rm(list=ls())
 at the start of each file and then set the seed to 42. It does not affect 
 any conclusions, but may be the reason if you find you get a slightly 
 different set of random numbers (and thus fitted values of b) to those
 in my provided .RData files.  I am running R in an Emacs shell through
 Emacs Speaks Statistics, which may or may not be important. Further details
 on this are in readRerun.txt, which is really just my ongoing notes to 
 document the issue.

