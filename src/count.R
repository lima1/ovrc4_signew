#load("input/eset.all.scaled.rda")
esets.all <- esets[-23]

esets.all$GSE19829.GPL8300_eset$histological_type <- "ser"
sum(sapply(esets.all, ncol))

XX <- lapply(esets.all, function(X) X[,X$histological_type=="ser"])
sum(sapply(XX, ncol))

XX <- sapply(XX, function(X) X[,is.na(X$summarystage) |
    X$summarystage=="late"])
XX <- sapply(XX, function(X) X[,is.na(X$summarygrade) |
    X$summarygrade=="high"])
XX$TCGA_eset_early <- esets.validation$TCGA_eset
length(XX)
sum(sapply(XX, ncol))
XX <- lapply(XX, function(X) 
X[,!is.na(X$days_to_death) | !is.na(X$os_binary)])

sum(sapply(XX, ncol))
length(XX)

XX <- XX[sapply(XX, ncol)>=40 |
         sapply(XX, function(X) sum(X$summarystage=="early",na.rm=TRUE))>0 |
         sapply(XX, function(X) sum(!is.na(X$os_binary))>0) 
]

sum(sapply(XX, ncol))
length(XX)
XX <- XX[sapply(XX, function(X) sum(X$vital_status=="deceased"))>=15 | 
         sapply(XX, function(X) sum(X$summarystage=="early",na.rm=TRUE))>0 |
         sapply(XX, function(X) sum(!is.na(X$os_binary))>0) 
]

XX$PMID15897565_eset <- NULL
sum(sapply(XX, ncol))

X <- XX$GSE17260_eset 
Y <- XX$GSE32062.GPL6480_eset
isc <- intersect(gsub("s","", X$alt_sample_name), gsub("d.*","",
Y$alt_sample_name))
iscx <- sapply(isc, function(x) grep(paste("^",x,"d",sep=""),
 Y$alt_sample_name))
XX$GSE32062.GPL6480_eset <- XX$GSE32062.GPL6480_eset[,-iscx]
XX <-.fixTCGA(XX)
sum(sapply(XX, ncol))
length(XX)

XXX <- c(esets.f, esets.validation, esets.binary[1])

sum(sapply(XXX, ncol))
