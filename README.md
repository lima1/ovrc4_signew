# Risk Prediction for late-stage ovarian cancer by meta-analysis of 1,525 patient samples: Biologic and Clinical Correlations

To reproduce the results of the
paper, simply fetch the entire repository using Mercurial with the
command below (or download from the "download" link at the url
<https://bitbucket.org/lima1/ovrc4_signew>):

	hg clone https://lima1@bitbucket.org/lima1/ovrc4_signew


Enter the newly created directory: 

     cd ovrc4_signew

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

    * **metaCMA.R**: code of our algorithm

    * **utils.R**: small helper functions, moved here to make the main
     analysis code more readable.

    * **ggplot2extras.R**: extends ggplot2 to make it look more like normal 
      R plots.   

    * **count.R**: counts the samples as shown in Figure 1 


* **input/** - input files

    * **runEsetList.sh**: create ExpressionSet objects from curatedOvarianData

    * **patientselection2...**: various patient filters for curatedOvarianData

    * **TCGA_ovsig.RData**: The TCGA model as survHD object 

    * **qrtpcr.xls**: The qRT-PCR data presented in the paper

    * **IHC.xls**: The IHC data presented in the paper

    * **bonometcgaduplicates.xls**: Spreadsheet with Bonome samples included
        in TCGA

    * **JCI65833sd1.xls**: The Verhaak et al. Supplemental Table with their
        signature    
        
    * **TCGA_489_UE.k4.txt**: The TCGA subtypes

    * **AOCS_subtytpes.xlsx**: Spreadsheet with the AOCS subtypes

    * **berchuck04.csv**: The manually extracted Berchuck 2004 signature.

------------------------------------------------------
------------------------------------------------------
# Requirements #

R/Bioconductor: Tested on R 3.0.1 and Bioconductor 2.12, but the
code should work with any relatively recent versions. The texreg package is
required in fairly recent version.

All other needed packages will be installed automatically when running
the install scripts shown above.

