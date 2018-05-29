readMeBinStructure.txt - description of R code to reproduce results from second
manuscript (on analysing binned data) and for analysing new such data, for
submission to ICES Journal of Marine Science.

Figure numbers are those upon submission. So numbers may need changing here at some
point. I'm taking my original Sweave files and just extracting the necessary
code to reproduce each figure (let me know if you would like the original
Sweave files, which automates some text also).

Contents and how to create each figure and results table
********************************************************

countsFunctions.r - common functions called by other code.

Length-weight relationships and consequences
--------------------------------------------

binFigs/binLW.r
Figure 1 - Length-weight relationships for two species.
Figure 2 - how binned body masses get assigned to logarithmic bins.


Simulated data
--------------

Figure 3 (histograms of b for simulated data sets) and Figure 4 (confidence
 intervals of b for simulated data sets) and Table A.3 (summary statistics
 related to those results).
binFitting/fitMLEmidMLEbin.r - create the simulated data
binFitting/fitMLEbinConf.r - does Figures 3, 4, Table A.3 and extra figures
 for each method.

binFitting/xmin16/ - repeating simulations but for xmin=16 for Figures A.36 and
 A.37 and Table A.4.
binFitting/xCutOff16/ - repeating simulations but for a cut-off of 16 for
 Figures A.38 and A.39 and Table A.5.


North Sea IBTS data
-------------------

nSeaFung/nSeaFungData/ibtsQ1cpuelength.RData - original file of data as downloaded
 by Julia Blanchard from ICES DATRAS website.

nSeaFung/nSeaFungData/speccodes.csv - species codes for the data.


nSeaFung/nSeaFungImport.Snw - Sweave file to load data, understand it and tidy
 up. Run using Sweave("nSeaFungImport.Snw") - this runs the R code and creates
 nSeaFungImport.tex latex file, which you can then run
 (e.g. pdflatex nSeaFungImport.tex)
 to create a .pdf. If you don't use latex then R code will still have run, done
 the tidying up, and created nSeaFungImport.Rdata.

nSeaFung/nSeaFungAnalysis.Snw - Sweave file to show example data (Table 1),
 fit all eight methods (from MEE paper) to each year of the IBTS data and fit
 linear regressions to the estimated size-spectrum exponents (Figure A.1 and Table
 A.1 [except MLEbins row - see below]). 

nSeaFung/nSeaMLEbins.Snw - Sweave file to analyse the IBTS data using the
 MLEbins method. Creates Figure 5 (and related Figures A.2, A.3 and A.4) showing
 species-specific body mass bins resulting from the length bins, Figure 7
 (comparison of MLE and MLEbins values of b through time) and MLEbins row of
 Table A.1. 

nSeaFung/nSeaMLEbins-recommend.Snw - recommended plots of data and fits, Figures
 6 and A.6-A.35, associate Table A.2 of results, plus Figure A.5
 (two horizontal red and pink bars) to show how we obtain the ranges of
 counts. [Adapting from original nSeaFungMLEbinsRecommend-ISD.Snw.]


Example histograms
------------------

exampleHists/histBinTypes.r - creates example histograms to show effects of
 binning, Figures A.40-A.44.

