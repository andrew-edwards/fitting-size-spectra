# fitting-size-spectra

This repository contains R code for the manuscript 

**Testing and recommending methods for fitting size spectra to data** 

by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard

that we submitted and are currently revising for *Methods in Ecology and Evolution*. 

I will update this information if the manuscript gets accepted. Meanwhile I am migrating my code for the manuscript from a private GitHub repository to this public one, so that I can improve the code and then just zip up a version to submit to the journal. I am also improving the documentation so so that others can use the code to analyse their own data and repeat and extend our analyses.

To download the code from the GitHub site just click the 'Download ZIP' button (near the top on the right). If you use GitHub then feel free to fork and even adapt the code. 

This repository is **under development** and will change quite a bit while I'm finalising the revisions, though really just in the file structure and documentation. The code all currently works (though it may not yet run properly because I've just changed the file structure). Let me know if you want to use the code before I finish the revisions (end of May 2016).

# Repository Contents

**README.md** - this file

**.gitignore** - ignore this if if you don't know what it is

**code/** - directory containing all the R code

**code/readMeCode.txt** - readme file for the code directory. Contains instructions and details of the files and further subdirectories.

The subdirectories of **code/** are summarised below, but see **readMeCode.txt** for full details.

**code/single/** - testing the eight methods on a single data set.

**code/multiple/** - testing the eight methods on multiple (10,000) simulated data sets. Contains subdirectories for the sensitivity tests.

**code/MLEbin/** - MLEbin method for likelihood when the data are only available in binned form.

**code/recommend/** - recommended likelihood calculations and resulting plots of data and fitted size spectrum (Figure 6).
 