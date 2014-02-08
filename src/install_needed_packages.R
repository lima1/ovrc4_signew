cran.packages <- c("xtable", "gplots", "TeachingDemos", "logging",
                   "survival", "rmeta", "survcomp", "snow", "rlecuyer",
                   "RColorBrewer", "HGNChelper", "metafor", "ROCR", 
                   "pROC", "maxstat", "lattice", "texreg",
                   "CoxBoost", "gbm", "survIDINRI",
                   "randomSurvivalForest", "uniCox", "superpc", "pensim",
                   "RUnit", "cvTools", "ez", "sampling", "impute", "knitr", "devtools")

bioc.packages <- c("genefilter", "affy", "preprocessCore", "hgu133a.db",
    "limma", "GSVA", "impute", "survcomp", "GEOquery", "hu6800.db",
    "multtest","graphite", "genomes"
    )

for (pkg in cran.packages){
    if (!require( pkg, character.only=TRUE ))
        install.packages(pkgs=pkg, repos="http://cran.us.r-project.org")
}

if (!require(BiocInstaller))
    stop("You need to install Bioconductor, which includes BiocInstaller.")

for (pkg in bioc.packages){
    if (!require( pkg, character.only=TRUE ))
        biocLite(pkgs=pkg, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
}

package.name <- "curatedOvarianData"   ##alternatives are NormalizerVcuratedOvarianData and FULLVcuratedOvarianData
package.version <- "0.99.21"

library(devtools)

##Do the installation:
package.url <- paste("http://bcb.dfci.harvard.edu/ovariancancer/dfiles_old/", package.name, "_",
                     package.version, ".tar.gz", sep="")
if( !require(package.name, character.only=TRUE) || package.version(package.name) != package.version ){
    install_url(package.url)
}

install_github("LeviRmisc", user="lwaldron")

