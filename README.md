# Reproducible analysis for the validation of published signatures of ovarian cancer survival #

To reproduce the results of the
paper, simply fetch the entire repository using Mercurial with the
command below (or download from the "download" link at the url
<https://bitbucket.org/lima1/ovrc4_signew>):

	hg clone ssh://hg@bitbucket.org/lima1/ovrc4_signew

Enter the newly created directory and start R: 

     cd ovrc4_sigvalidation
     R
     
Run the analysis with knitr in R:
    
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

    * **createEsets.sh**: create ExpressionSet objects from curatedOvarianData

    * **TCGA_ovsig.RData**: The TCGA model as survHD object 

------------------------------------------------------
------------------------------------------------------
# Requirements #

R/Bioconductor: Tested on R 2.15.0 and Bioconductor 2.10, but the
pipeline should work on any relatively recent versions.

All other needed packages will be installed automatically when running
the pipeline, by the following scripts in src/:

* install_needed_packages.R
* install_survHD.R
* createEsets.R   (installs curatedOvarianData)

