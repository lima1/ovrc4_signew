# Risk prediction by meta-analysis for late stage ovarian cancer

To reproduce the results of the
paper, simply fetch the entire repository using Mercurial with the
command below (or download from the "download" link at the url
<https://bitbucket.org/lima1/ovrc4_signew>):

	hg clone ssh://hg@bitbucket.org/lima1/ovrc4_signew


Enter the newly created directory: 

     cd ovrc4_sigvalidation

Install dependencies in R:

    R    
    source("src/install_needed_packages.R")
    source("src/install_survHD.R")
    q()

Fetch the data from the curatedOvarianData package:

    cd input
    ./runEsetList.sh
    cd ..

Run the analysis with knitr in R:

    R    
    library(knitr)
    knit("metasig.Rnw")

------------------------------------------------------
------------------------------------------------------
# Summary of files and directories

* **src/**  - additional source code

    * **utils.R**: small helper functions, moved here to make the main
     analysis code more readable.

    * **metaCMA.R**: code of our algorithm

* **input/** - input files

    * **runEsetList.sh**: create ExpressionSet objects from curatedOvarianData

    * **patientselection2...**: various patient filters for curatedOvarianData

    * **TCGA_ovsig.RData**: The TCGA model as survHD object 

------------------------------------------------------
------------------------------------------------------
# Requirements #

R/Bioconductor: Tested on R 2.15.1 and Bioconductor 2.11, but the
pipeline should work on any relatively recent versions.

All other needed packages will be installed automatically when running
the install scripts shown above.

