cran.packages <- c("xtable", "gplots", "TeachingDemos", "logging",
                   "survival", "rmeta", "survcomp", "snow", "rlecuyer",
                   "RColorBrewer", "HGNChelper", "metafor", "ROCR" )

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
