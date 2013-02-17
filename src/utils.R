# datasets we use for training the debulking prediction model
.debulkingFilter <- function(X) { 
     ncol(X) < 50 | sum(X$debulking=="suboptimal",na.rm=TRUE) == 0 
}

# takes an ExpressionSet from curatedOvarianData and formats a dataset name
.getDatasetNames <- function(esets) {
    s <- sapply(esets, function(Y) gsub("Cancer GenomeAtlas Research Network", 
        "TCGA", experimentData(Y)@lab))
    s <- gsub(",.*20", " et al, 20", s)
    gsub(" hgu95", "", s)
}

# adds a label to each ExpressionSet, classifying each sample as long or
# short-term survivor
.dichotomizeshortlong <- function(eset, s=365.25, l=365.2*4,
label="os_my_binary") {
    eset <- eset[, eset$y[,1] > s  |
        eset$y[,2] == 1]

    eset[[label]] <- as.factor(ifelse(eset$y[,1] < s, "short",
    ifelse(eset$y[,1] > l,"long", NA)))
    eset[, !is.na(eset[[label]])]
}

.dichotomizeshortlongdf <- function(eset, s=365.25, l=365.2,
label="os_1yr") {
    eset <- eset[eset$y[,1] > s  |
        eset$y[,2] == 1, ]

    eset[[label]] <- as.factor(ifelse(eset$y[,1] < s, "short",
    ifelse(eset$y[,1] > l,"long", NA)))
    eset[!is.na(eset[[label]]),]
}

# Display the dependency of a parameter on the prediction accuracy yi
.plotN <- function(esets, yi, xlab="Number of Genes",
ylab="C-Index (Concordance)") {
    w <- sapply(esets,ncol)
    dfCV <- stack(as.data.frame(do.call(rbind, yi)))
    dfCV <- cbind(dfCV, Genes=c(5,10,seq(25,250,25)))
    dfCV$ind <- unlist(lapply(.getDatasetNames(esets), rep, 12))
    ggplot(dfCV, aes(Genes, values))+geom_line()+facet_wrap(~ind) +
    xlab(xlab) + 
    ylab(ylab) + 
    theme_classic2(base_size=13) + 
    theme(axis.text.x=element_text(angle=45,hjust=1))
}

# show a Kaplan-Meier analysis of a leave-one-dataset-out cross-validation
.lodocvPlot <- function(X, models, ids = 1:length(X), nrow = 2, ncol = 3,
censor.at = 365.25 * 5, cutpoint=NULL, plot=TRUE,...) {
#    X <- X[-match("TCGA_eset", names(X))]
    i <- !sapply(X, .defaultFilter)
    names(X)[i] <- .getDatasetNames(X[i])
    names(X)[i] <- paste(1:length(X[i]),". ", names(X)[i], sep="")
    par(mfrow = c(nrow, ncol))
    par(mar = c(4.3, 4.5, 3, 1))
    res <- lapply(ids, function(i) {
        # calculate median risk score in the training data
        if (length(models) == 1) model <- models
        else model <- models[[i]]
        if (is.null(cutpoint))
            cutpoint <- median(unlist(lapply(X[-i], function(x) predict(model, newdata=x)@lp)))
        local.plot <- plot 
        if (.defaultFilter(X[[i]])) local.plot=FALSE
        plotKMStratifyBy(cutpoints = cutpoint,
        linearriskscore = predict(model,newdata=X[[i]])@lp, y = X[[i]]$y, 
        censor.at = censor.at, cex.base = 1.4, show.n.risk = FALSE, show.HR = FALSE, 
        show.legend = FALSE, main = names(X)[i], plot = local.plot, ...)
    })
    for (i in 1:length(res)) res[[i]]$y <- X[ids][[i]]$y
    res
}

.valPlot <- function(model, coefficients = names(model@coefficients), plot = TRUE) {
    # use the median of the training risk scores for the Japanese and TCGA data,
    # for the other the median of the test data because their platforms have not
    # all genes of the signature
    risk.fmtrain <- lapply(esets.f, function(X) predict(final.model, newdata = X)@lp)
    cutoff <- median(unlist(risk.fmtrain))
    esets.f.too.small <- esets.f[sapply(esets.f, .defaultFilter)]
    
    titles <- c(
        paste(c("A) ", "B) "), .getDatasetNames(esets.f.too.small)),
        "C) Konstantinopoulos 2010",
        paste("D) Gillet 2012 (", sum(coefficients %in% featureNames(esets.validation$GSE30009_eset)), 
        " genes)",sep=""), 
        paste("E) ",
            .getDatasetNames(list(esets.validation$GSE32062.GPL6480_eset))), 
        "F) Early Stage TCGA")
    cutpoints <- list(cutoff, cutoff, NULL, NULL, cutoff, cutoff)
    
    par(mfrow = c(3, 3))
    par(mar=c(4.5, 4.1, 2.5, 1.5))
    res <- lapply(1:length(esets.f.too.small), function(i) plot(model, newdata = esets.f.too.small[[i]], 
        newy = esets.f.too.small[[i]]$y, show.n.risk = FALSE, show.legend = FALSE, 
        show.HR = FALSE, cex.base = 1.4, main = titles[i], cutpoints = cutpoints[[i]], 
        plot = plot, censor.at = 365.25 * 5))
    res <- lapply(1:length(esets.validation), function(i) plot(model, newdata = esets.validation[[i]], 
        newy = esets.validation[[i]]$y, show.n.risk = FALSE, show.legend = FALSE, 
        show.HR = FALSE, cex.base = 1.4, main = titles[i+2], cutpoints =
        cutpoints[[i+2]], 
        plot = plot, censor.at = 365.25 * 5))
    for (i in 1:length(esets.validation)) res[[i]]$y <- esets.validation[[i]]$y
    
    # Berchuck is a subset of Dressman 
    for (i in 1:1) {
        pred <- predict(model, newdata = esets.binary[[i]], type = "lp")@lp
        paneln <- c("G)", "H)")
        if (plot) 
            res.roc <- .plotROC(pred, esets.binary[[i]]$os_binary, main
            = paste(paneln[i],.getDatasetNames(esets.binary)[i]))
    }
    res
}

.lodocvPlotForest <- function(X, main = "C-Index") {
    if ("TCGA_eset" %in% names(X)) 
        X <- X[-match("TCGA_eset", names(X))]

    names(X) <- .getDatasetNames(X)

    metaCMA.forest.label(esets = X, label = "risk", y = "y", concordance = TRUE, 
        cex = 0.7, main = main, xlab = "")
}


.plotROC <- function(pred, labels, plot = TRUE, na.rm = TRUE, colorize = FALSE, ...) {
    require(ROCR)
    require(pROC)
    if (na.rm) {
        idx <- !is.na(labels)
        pred <- pred[idx]
        labels <- labels[idx]
    }
    pred.rocr <- prediction(pred, labels)
    perf.rocr <- performance(pred.rocr, "tpr", "fpr")
    auc <- performance(pred.rocr, "auc")@y.values[[1]][[1]]
    roc.obj <- roc(labels, pred)
    auc.ci <- ci(roc.obj)
    significant <- ifelse(ci(roc.obj, conf.level=0.9)[1] > 0.5, "*", "")
    best <- coords(roc.obj,x="best") 
    if (plot) {
        plot(perf.rocr, colorize = colorize, cex.lab = 1.3, ...)
        text(0, 0.9, paste("AUC = ", round(auc, digits = 2), significant,
        sep=""), cex = 1.5, pos = 4)
        abline(a = 0, b = 1, lty = 2)
        text(1, 0.1, paste("n =", length(labels)), cex = 1.5, pos = 2)
        abline(a = 0, b = 1, lty = 2)
    }
    invisible(list(auc,auc.ci,best))
}


.plotROCpanel <- function(preds, labels, titles, nrow, ncol,...) {
    par(mfrow = c(nrow, ncol))
    aucs <- lapply(1:length(preds), function(i)
        .plotROC(preds[[i]], labels[[i]], main = titles[[i]],...))

    # calculate the pooled AUC     
    sei <- sapply(aucs, function(auc) (auc[[2]][2]-auc[[2]][1])/(3.92/2))
    yi <- sapply(aucs, function(auc) (auc[[1]]))
    list(rma=metafor::rma(yi=yi,sei=sei,method="FE"),optimal.cutpoints=sapply(aucs,
    function(auc) auc[[3]]),yi=yi,sei=sei)
}

.forestPlotDebulking <- function(preds1, labels, prior=NULL, titles, method="FE", ...) {
    dat <- data.frame(opt=sapply(labels, function(x) sum(x=="optimal")),
    subopt=sapply(labels, function(x) sum(x=="suboptimal")))
    
    # get optimal cutpoints
    best <-  .optimalROCclassification(preds1, labels, best.method="closest.topleft")[[1]]
    
    # classify samples according cutpoints (either best or the ones keeping the
    # correct ratio in the validation data) 
    classification <- lapply(1:length(preds1), function(i) {
        p1scaled <- scale(preds1[[i]])
        if (!is.null(prior)) p1scaled <- as.factor(p1scaled > quantile(p1scaled,
            p=1-prior[i]))
        else  p1scaled <- as.factor(p1scaled > best[i])
        p1scaled
    })

    res <- lapply(1:length(preds1), function(i) {

    summary(glm(labels[[i]]~classification[[i]],
    family="binomial"))$coefficients})
    #res1 <- metafor::rma(yi, vi, data=dat, method=method)
    yi <- sapply(res, function(x) x[2,1])
    sei <- sapply(res, function(x) x[2,2])
    rma1 <- metafor::rma(yi=yi, sei=sei, method=method)

    forest.rma(rma1,
    atransf=exp,xlim=c(-3,4),ilab=dat,
    at= log(sapply(-2:5, function(x) 1*2.1^x)),
    ilab.xpos=c(-6.5,-4.5)*3.5/16,slab=titles,mlab="Overall",
    xlab="Odds Ratio (log scale)")
    op <- par(font=2)
    text( -3, 10, "Dataset",pos=4)
    text(c(-6.5,-4.5)*3.5/16,10,c("Opt.", "Subopt."))
    text(c(-5.5)*3.5/16,11,c("Debulking"))
    text(4,10, "Gene Signature Odds Ratio [95% CI]",pos=2)
    par(op)

    list(rma1,classification,dat)
}

.reclassPlot <- function(y, r1, r2, t0, title, plot = TRUE) {
    x <- IDI.INF(y, r1, r2, t0 = t0, npert = 600)
    y <- IDI.INF.OUT(x)
    if (plot) {
        IDI.INF.GRAPH(x, main = paste(title, " (t0 = ", t0, ")", sep = ""), las = 1, 
            cex.lab = 1.3, cex.axis=1.3)
        text(par("usr")[1], 0.95, paste("IDI = ", y[1, 1], "; 95% CI ", y[1, 2], 
            " to ", y[1, 3], sep = ""), pos = 4, cex = 1.45)
    }
    y
}

# create a data.frame with clinical covariates out of a list of ExpressionSets
# from curatedOvarianData
.createClinical <- function(esets, clinical.covars = NULL) {
    if (is.null(clinical.covars)) {
        clinical.covars <- varLabels(esets[[1]])
        clinical.covars <- clinical.covars[clinical.covars != "y"]
    }        
    require(plyr)
    df <- ldply(esets, function(X) as(phenoData(X),
        "data.frame")[,clinical.covars])

    df$batch <- as.factor(unlist(sapply(esets, function(X) rep(experimentData(X)@lab, 
        ncol(X)))))
    
    df$debulking <- as.factor(gsub("unknown", NA, df$debulking))
    df$y <- Surv(as.numeric(df$days_to_death), df$vital_status == "deceased")
    df
}

# calculate the p-value of overlap for two gene signatures
.pOverlap <- function(A, B, background) {
    n1 <- length(A)
    n2 <- length(B)
    n <- length(background)
    m <- length(intersect(A, B))
    phyper(min(n1, n2), n1, n - n1, n2) - phyper(m - 1, n1, n - n1, n2)
}

# given a model that uses gene symbols, translate it to a model using
# probesets. works only with curatedOvarianData ExpressionSets, because we
# assume a maxmean_probeset slot here
.modelSymbolToAffy <- function(model, ref.eset) {
    model.affy <- model
    names(model.affy@coefficients) <-
    featureData(ref.eset[names(model@coefficients),])$maxmean_probeset
    model.affy
}

.boostrapHRs <- function(data, r1="strata", r2="strata_tcga", n=500,
inverse=FALSE) {
    .doBS <- function(i) {
        idx <- sample(nrow(data), replace=TRUE)
        ldata <- data.frame(y=Surv(data[,1], data[,2]), r1=data[[r1]], r2=data[[r2]])
        ldata <- ldata[idx,]
        colnames(ldata) <- c("y", "r1", "r2")

        f1 <- coxph(y~r1,ldata)
        f2 <- coxph(y~r2,ldata)
        if (inverse) hrdiff <- 1/summary(f1)$conf.int[1] - 1/summary(f2)$conf.int[1]
        else hrdiff <- summary(f1)$conf.int[1] - summary(f2)$conf.int[1]
        hrdiff
    }
    sapply(1:n, .doBS)
}


# finds gene sets predictive of the specified label,
# simply uses limma for differential expression
.getFingerprintGeneSets <- function(eset, labels,
p.value,lfc,number) {
    require(limma)
    .calcLimma <- function(l) {
        TS <- rep("Ctrl", length(labels))
        TS[labels==l] <- "Test"
        TS <- as.factor(TS)
        design <- model.matrix(~0+TS)
        colnames(design) = levels(TS)
        fit <- lmFit(exprs(eset), design = design)
        cont <-
            makeContrasts(Test-Ctrl,levels=design)
        fit2 <- contrasts.fit(fit, cont)
        fit2 <- eBayes(fit2)
    }
    limma.res <- lapply(levels(labels), .calcLimma)
    filtered.res <- lapply(limma.res, function(x) topTable(x,
    number=5000,p.value=p.value,lfc=lfc)[,1:2])
    genesets <- lapply(filtered.res, function(x)
    list(Up=utils::head(x[x[,2]>0,1],number),
        Down=utils::head(x[x[,2]<0,1],number)))
    for (i in 1:length(genesets)) {
        names(genesets[[i]]) <- paste(levels(labels)[i], names(genesets[[i]]),
        sep="_")
    }    
    genesets <- unlist(genesets, recursive=FALSE)
    genesets[order(gsub(".*_","",names(genesets)),decreasing=TRUE)]
}

.classifyGSVA <- function(X=NULL, eset, genesets,...) {
    if (is.null(X)) X <- gsva(exprs(eset),
    genesets,mx.diff=TRUE,parallel.sz=1,...)[[1]]
    Y <- X[1:4,] - X[5:8,]
    rownames(Y) <- gsub("_Up","", rownames(Y))
    Y
}

.optimalROCclassification <- function(preds, labels, ...) {
    .doIt <- function(i) {
        roc.obj <- roc(unlist(labels[-i]), unlist(preds[-i]))
        coords(roc.obj, x="best", ...)[1]
    }
    best <- sapply(1:length(preds), .doIt)
    classification  <- lapply(1:length(best), function(i) preds[[i]] > best[i])
    cohorts <- lapply(1:length(preds), function(i) table(classification[[i]], labels[[i]]))
    combined <- table(unlist(classification), unlist(labels))
    list(best, cohorts, combined)
}

cvRisk <- function(fit, data, y=data$y, linear=TRUE,...) {
    folds <- cvFolds(nrow(data),...)
    yhat <- cvFit(fit, y=y, data=data, folds=folds, cost=function(y, yhat) yhat )
    folds.order <- folds$subset[order(folds$which)]
    if (linear) return(yhat$cv[order(folds.order)]  )
    list(risk=lapply(1:max(folds$which), function(i)
    yhat$cv[sort(folds$which)==i]),
         y   =  lapply(1:max(folds$which), function(i)
         y[folds$subset[folds$which==i] ]))
}

