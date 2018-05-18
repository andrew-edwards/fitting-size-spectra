readMeBinStructure.txt - description of R code to reproduce results from second
manusript (on analysing binned data) and for analysing new such data.

Figure numbers are those from version ssmPart2-2018-04-13.pdf, which hopefully
shouldn't change before submission. So numbers may need changing here at some
point. I'm taking my original Sweave files and just extracting the necessary
code to reproduce each figure (let me know if you would like the original
Sweave files, which automates some text also).


countsFunctions.r - common functions called by other code.

Figure 1 - Length-weight relationships for two species.
binFigs/binLW.r

Figure 2 - how binned body masses get assigned to logarithmic bins.
binFigs/binLW.r

Figure 3 (histograms of b for simulated data sets) and Figure 4 (confidence
 intervals of b for simulated data sets) and Table A.3 (summary statistics
 related to those results).
binFitting/fitMLEmidMLEbin.r - create the simulated data
binFitting/fitMLEbinConf.r - does Figures 3, 4, Table A.3 and extra figures
 for each method.

binFitting/xmin16/ - repeating simulations but for xmin=16

North Sea IBTS data:

nSeaFung/nSeaFungData/ibtsQ1cpuelength.RData - original file of data as downloaded
 by Julia Blanchard from ICES DATRAS website.

nSeaFung/nSeaFungData/speccodes.csv - species codes for the data.


nSeaFung/nSeaFungImport.Snw - Sweave file to load data, understand it and tidy
 up. Run using Sweave("nSeaFungImport.Snw") - this runs the R code and creates
 nSeaFungImport.tex latex file, which you can then run
 (e.g. pdflatex nSeaFungImport.tex)
 to create a .pdf. If you don't use latex then R code will still have run, done
 the tidying up, and created nSeaFungImport.Rdata.

nSeaFung/nSeaFungAnalysis.Snw - Sweave file to show example data, fit all eight
 methods (from MEE paper) to each year of the IBTS data, and fit linear
 regressions to the estimated size-spectrum exponents. 

nSeaFung/nSeaMLEbins.Snw - Sweave file to analyse the IBTS data using the
 MLEbins method. Also creates figures...   .
 [Needs some of nSeaFungMLEbin.Snw (because that does some pre-analysis),
 nSeaFungMLEbinsNew.Snw, and nSeaFungCompareNew.Snw].

nSeaFung/nSeaMLEbins-recommend.Snw - recommended plots plus the schematic
 diagram (two horizontal red and pink bars) to show how we obtain the ranges of
  counts. Adapting from original nSeaFungMLEbinsRecommend-ISD.Snw.

Figure 5 - species-specific body mass bins resulting from the length bins


