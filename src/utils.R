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

# show a Kaplan-Meier analysis of a leave-one-dataset-out cross-validation
.lodocvPlot <- function(X, models, ids = 1:length(X), nrow = 3, ncol = 3, censor.at = 365.25 * 5, ...) {
#    X <- X[-match("TCGA_eset", names(X))]
    names(X) <- .getDatasetNames(X)
    names(X) <- paste(1:length(X),". ", names(X), sep="")
    par(mfrow = c(nrow, ncol))
    par(mar = c(4.3, 4.5, 3, 1))
    res <- lapply(ids, function(i) {
        # calculate median risk score in the training data
        if (length(models) == 1) model <- models
        else model <- models[[i]]
        cutpoint <- median(unlist(lapply(X[-i], function(x) predict(model,
        newdata=x)@lp)))

        plotKMStratifyBy(cutpoints = cutpoint,
        linearriskscore = predict(model,newdata=X[[i]])@lp, y = X[[i]]$y, 
        censor.at = censor.at, cex.base = 1.4, show.n.risk = FALSE, show.HR = FALSE, 
        show.legend = FALSE, main = names(X)[i], ...)
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
    
    titles <- c(
        "A) Konstantinopoulos 2010",
        paste("B) Gillet 2012 (", sum(coefficients %in% featureNames(esets.validation$GSE30009_eset)), 
        " genes)",sep=""), 
        paste("C) ",
            .getDatasetNames(list(esets.validation$GSE32062.GPL6480_eset))), 
        "D) Early Stage TCGA")
    cutpoints <- list(NULL, NULL, cutoff, cutoff)
    
    par(mfrow = c(2, 3))
    res <- lapply(1:length(esets.validation), function(i) plot(model, newdata = esets.validation[[i]], 
        newy = esets.validation[[i]]$y, show.n.risk = FALSE, show.legend = FALSE, 
        show.HR = FALSE, cex.base = 1.4, main = titles[i], cutpoints = cutpoints[[i]], 
        plot = plot, censor.at = 365.25 * 5))
    for (i in 1:length(esets.validation)) res[[i]]$y <- esets.validation[[i]]$y
    
    
    for (i in 1:length(esets.binary)) {
        pred <- predict(model, newdata = esets.binary[[i]], type = "lp")@lp
        if (i == 1) {
            Xtmp <- lapply(esets.f, .dichotomizeshortlong, 365.25 * 3, 365.25 * 7)
        } else {
            Xtmp <- lapply(esets.f, .dichotomizeshortlong, 365.25 * 5, 365.25 * 5)
        }
        pred.logit <- .p2logitSingle(Xtmp, pred, model, y = "os_my_binary")
        paneln <- c("E)", "F)")
        if (plot) 
            res.roc <- .plotROC(pred.logit, esets.binary[[i]]$os_binary, main
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


.plotROC <- function(pred, labels, plot = TRUE, na.rm = TRUE, ...) {
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
    auc.ci <- ci(roc(labels,pred))
    
    if (plot) {
        plot(perf.rocr, colorize = TRUE, cex.lab = 1.3, ...)
        text(0, 0.9, paste("AUC =", round(auc, digits = 2)), cex = 1.5, pos = 4)
        abline(a = 0, b = 1, lty = 2)
        text(1, 0.1, paste("n =", length(labels)), cex = 1.5, pos = 2)
        abline(a = 0, b = 1, lty = 2)
    }
    invisible(list(auc,auc.ci))
}

.plotROCpanel <- function(preds, labels, titles, nrow, ncol,...) {
    par(mfrow = c(nrow, ncol))
    aucs <- lapply(1:length(preds), function(i)
        .plotROC(preds[[i]], labels[[i]], main = titles[[i]],...))

    # calculate the pooled AUC     
    sei <- sapply(aucs, function(auc) (auc[[2]][2]-auc[[2]][1])/(3.92/2))
    yi <- sapply(aucs, function(auc) (auc[[1]]))
    metafor::rma(yi=yi,sei=sei,method="FE")
}

.p2logitSingle <- function(esets.train, preds, model, y = "os_my_binary") {
    esets.c <- .combineEsets(esets.train, y = y)
    tmp.df <- data.frame(risk = predict(model, newdata = esets.c$X)@lp, covar = esets.c$y)
    
    fit <- glm(covar ~ risk, tmp.df, family = "binomial")
    predict(fit, newdata = data.frame(risk = preds), type = "response")
}

.p2logit <- function(esets, fits, preds = NULL, y = "debulking") {
    .doIt <- function(i) {
        esets.c <- .combineEsets(esets[-i], y = y)
        tmp.df <- data.frame(risk = predict(fits[[i]]$fit, newdata = esets.c$X)@lp, 
            covar = esets.c$y)
        fit <- glm(covar ~ risk, tmp.df, family = "binomial")
        risk <- fits[[i]]$risk@lp
        if (!is.null(preds)) 
            risk <- preds[[i]]
        predict(fit, newdata = data.frame(risk = risk), type = "response")
    }
    lapply(1:length(esets), .doIt)
}

.p2logitMV <- function(esets, fits, preds = NULL, y = "debulking", x1 = "tumorstage") {
    .doIt <- function(i) {
        esets.c <- .combineEsets(esets[-i], y = y)
        tmp.df <- data.frame(risk = predict(fits[[i]]$fit, newdata = esets.c$X)@lp, 
            covar = esets.c$y, x1 = unlist(sapply(esets[-i], function(X) X[[x1]])))
        fit <- glm(covar ~ risk + x1, tmp.df, family = "binomial")
        risk <- fits[[i]]$risk@lp
        if (!is.null(preds)) 
            risk <- preds[[i]]
        predict(fit, newdata = data.frame(risk = risk, x1 = esets[[i]][[x1]]), type = "response")
    }
    lapply(1:length(esets), .doIt)
}

.forestPlotDebulking <- function(preds1, labels, prior, titles, method="FE", ...) {
    dat <- data.frame(opt=sapply(labels, function(x) sum(x=="optimal")),
    subopt=sapply(labels, function(x) sum(x=="suboptimal")))

    res <- lapply(1:length(preds1), function(i) {
    p1scaled <- scale(preds1[[i]])
    summary(glm(labels[[i]]~p1scaled,
    family="binomial"))$coefficients})
    #res1 <- metafor::rma(yi, vi, data=dat, method=method)
    yi <- sapply(res, function(x) x[2,1])
    sei <- sapply(res, function(x) x[2,2])
    rma1 <- metafor::rma(yi=yi, sei=sei, method=method)

    forest.rma(rma1,
    atransf=exp,xlim=c(-3,2.5),ilab=dat,
    at= log(sapply(-2:5, function(x) 1*1.3^x)),
    ilab.xpos=c(-6.5,-4.5)*3.5/16,slab=titles,mlab="Overall",
    xlab="Odds Ratio (log scale)")
    op <- par(font=2)
    text( -3, 10, "Authors and Year",pos=4)
    text(c(-6.5,-4.5)*3.5/16,10,c("Opt.", "Subopt."))
    text(c(-5.5)*3.5/16,11,c("Debulking"))
    text(2.5,10, "Gene Signature Odds Ratio [95% CI]",pos=2)
    par(op)

    list(rma1,dat)
}

.reclassPlot <- function(y, r1, r2, t0, title, plot = TRUE) {
    x <- IDI.INF(y, r1, r2, t0 = t0, npert = 600)
    y <- IDI.INF.OUT(x)
    if (plot) {
        IDI.INF.GRAPH(x, main = paste(title, " (t0 = ", t0, ")", sep = ""), las = 1, 
            cex.lab = 1.3)
        text(par("usr")[1], 0.95, paste("IDI = ", y[1, 1], "; 95% CI ", y[1, 2], 
            " to ", y[1, 3], sep = ""), pos = 4, cex = 1.45)
    }
    y
}

.concPlotSize <- function(esets, ma, skip = 0, dataset_ids, ylab = "C-Index (Concordance)", 
    panel.num = "A)", ...) {
    labels <- .getDatasetNames(esets[dataset_ids])
    labels <- paste(labels, " (N = ", sapply(dataset_ids,function(i)
    ncol(esets[[i]])), ")",sep="")    
    # exclude TCGA, because we compare concordance to TCGA model throughout the
    # paper
    tcga_id <- match("TCGA_eset", names(esets))
    ids <- (1 + skip):length(esets)
    x <- sapply(ids, function(i) metaCMA.concordance(esets.f[-tcga_id], risks = lapply(ma[[i]]$fits[-tcga_id], 
        function(x) x$risk@lp))[[1]])
    
    lb <- unlist(x[5, ])
    ub <- unlist(x[6, ])
    
    par(mar = c(15, 5.5, 2, 0.5))
    plot(ids, unlist(x[1, ]), type = "l", las = 1, ylim = c(0.56, 0.65), xaxt = "n", 
        ylab = "", xlab = "", cex.axis = 1.1, cex.lab = 1.2, ...)
    axis(side = 1, at = ids, labels = labels[ids], las = 2, cex.axis = 1.2)
    title(ylab = ylab, cex.lab = 1.2, line = 4)
    points(ids, lb, col = "grey", lty = 2, type = "l")
    points(ids, ub, col = "grey", lty = 2, type = "l")
    mtext(panel.num, side = 3, cex = 1.5, at = 0)
    
}

# create a data.frame with clinical covariates out of a list of ExpressionSets
# from curatedOvarianData
.createClinical <- function(esets, clinical.covars = NULL) {
    if (is.null(clinical.covars)) 
        clinical.covars <- c("age_at_initial_pathologic_diagnosis", "debulking", 
            "tumorstage", "days_to_death", "vital_status")
    df <- data.frame(do.call(cbind, lapply(clinical.covars, function(covar) unlist(sapply(esets, 
        function(X) X[[covar]])))), stringsAsFactors = FALSE)
    colnames(df) <- clinical.covars
    df$tumorstage <- as.numeric(df$tumorstage)
    df$batch <- as.factor(unlist(sapply(esets, function(X) rep(experimentData(X)@lab, 
        ncol(X)))))
    
    df$debulking <- as.factor(gsub("unknown", NA, df$debulking))
    df$age_at_initial_pathologic_diagnosis <- as.numeric(df$age_at_initial_pathologic_diagnosis)
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

