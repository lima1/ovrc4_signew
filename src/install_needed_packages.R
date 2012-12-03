cran.packages <- c("xtable", "gplots", "TeachingDemos", "logging",
                   "survival", "rmeta", "survcomp", "snow", "rlecuyer",
                   "RColorBrewer", "HGNChelper", "metafor", "ROCR", 
                   "pROC", "maxstat", "lattice", "survIDINRI", "rockchalk", 
                   "cvTools")

bioc.packages <- c("genefilter", "affy", "preprocessCore", "hgu133a.db")

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
package.version <- "0.99.10"

##Do the installation:
package.url <- paste("http://bcb.dfci.harvard.edu/ovariancancer/dfiles_old/", package.name, "_",
                     package.version, ".tar.gz", sep="")
if( !require(package.name, character.only=TRUE) || package.version(package.name) != package.version ){
    library(devtools)
    install_url(package.url)
}

