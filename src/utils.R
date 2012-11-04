.lodocvPlot <- function(X, nrow = 3, ncol = 3, censor.at = 365.25 * 5, ...) {
    X <- X[-match("TCGA_eset", names(X))]
    
    names(X) <- sapply(X, function(Y) gsub("Cancer GenomeAtlas Research Network", 
        "TCGA", experimentData(Y)@lab))
    par(mfrow = c(nrow, ncol))
    par(mar = c(4.3, 4.5, 3, 1))
    res <- lapply(1:length(X), function(i) plotKMStratifyBy(cutpoints = median(unlist(sapply(X[-i], 
        function(x) x$risk)), na.rm = TRUE), linearriskscore = X[[i]]$risk, y = X[[i]]$y, 
        censor.at = censor.at, cex.base = 1.4, show.n.risk = FALSE, show.HR = FALSE, 
        show.legend = FALSE, main = names(X)[i], ...))
    for (i in 1:length(res)) res[[i]]$y <- X[[i]]$y
    res
}

.valPlot <- function(model, coefficients = names(model@coefficients), plot = TRUE) {
    # use the median of the training risk scores for the Japanese and TCGA data,
    # for the other the median of the test data because their platforms have not
    # all genes of the signature
    risk.fmtrain <- lapply(esets.f, function(X) predict(final.model, newdata = X)@lp)
    cutoff <- median(unlist(risk.fmtrain))
    
    titles <- c(paste("Gillet, 2012 (RT-PCR", sum(coefficients %in% featureNames(esets.validation$GSE30009_eset)), 
        "genes)"), "Konstantinopoulos 2010", experimentData(esets.validation$GSE32062.GPL6480_eset)@lab, 
        "Early Stage TCGA")
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
        
        if (plot) 
            res.roc <- .plotROC(pred.logit, esets.binary[[i]]$os_binary, main = experimentData(esets.binary[[i]])@lab)
    }
    res
}


.lodocvPlotForest <- function(X, main = "C-Index") {
    if ("TCGA_eset" %in% names(X)) 
        X <- X[-match("TCGA_eset", names(X))]
    
    names(X) <- sapply(X, function(Y) gsub("Cancer GenomeAtlas Research Network", 
        "TCGA", experimentData(Y)@lab))
    names(X) <- gsub(",.*", "", names(X))
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
    invisible(auc.ci)
}

.plotROCpanel <- function(preds, labels, titles, nrow, ncol) {
    par(mfrow = c(nrow, ncol))
    for (i in 1:length(preds)) {
        .plotROC(preds[[i]], labels[[i]], main = titles[[i]])
    }
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


.doTestAll <- function(preds, labels, priors, titles, method = "FE") {
    ci <- do.call(rbind, lapply(1:length(preds), function(i) try(.doTest(preds[[i]], 
        labels[[i]], prior = priors[i]))))
    rownames(ci) <- titles
    se <- do.call(rbind, lapply(1:length(preds), function(i) try(.doTestSE(preds[[i]], 
        labels[[i]], prior = priors[i]))))
    rownames(se) <- titles
    
    rma.res <- lapply(seq(1, 9, 2), function(i) metafor::rma(yi = se[, i], sei = se[, 
        i + 1], method = method))
    .fp <- function(p) round(p * 100, digits = 1)
    .fnp <- function(p) round(p, digits = 2)
    pooled1 <- lapply(rma.res[-5], function(x) c(.fp(x$b), paste(.fp(x$ci.lb), .fp(x$ci.ub), 
        sep = "-")))
    pooled2 <- lapply(rma.res[5], function(x) c(.fnp(x$b), paste(.fnp(x$ci.lb), .fnp(x$ci.ub), 
        sep = "-")))
    ci <- rbind(ci, unlist(c(pooled1, pooled2)))
}

.doTestSE <- function(preds, labels, prior, cutoff = quantile(preds, p = (1 - prior))) {
    idx <- !is.na(preds) & !is.na(labels)
    preds <- preds[idx]
    labels <- labels[idx]
    tbl <- table(preds > cutoff, labels)
    f <- fisher.test(as.factor(preds > cutoff), labels)
    .fp <- function(p) p
    .fnp <- function(p) p
    data.frame(PPV = .fp(tbl[2, 2]/sum(tbl[2, ])), PPV.SE = diff(binom.test(tbl[2, 
        2], sum(tbl[2, ]))$conf.int)/3.92, NPV = .fp(tbl[1, 1]/sum(tbl[1, ])), NPV.SE = diff(binom.test(tbl[1, 
        1], sum(tbl[1, ]))$conf.int)/3.92, FPR = .fp(tbl[2, 1]/sum(tbl[, 1])), FPR.SE = diff(binom.test(tbl[2, 
        1], sum(tbl[, 1]))$conf.int)/3.92, FNR = .fp(tbl[1, 2]/sum(tbl[, 2])), FNR.SE = diff(binom.test(tbl[1, 
        2], sum(tbl[, 2]))$conf.int)/3.92, OR = .fnp(f$estimate), OR.SE = diff(f$conf.int)/3.92, 
        stringsAsFactors = FALSE)
}

.doTest <- function(preds, labels, prior, cutoff = quantile(preds, p = (1 - prior))) {
    
    idx <- !is.na(preds) & !is.na(labels)
    
    preds <- preds[idx]
    labels <- labels[idx]
    tbl <- table(preds > cutoff, labels)
    print(tbl)
    f <- fisher.test(as.factor(preds > cutoff), labels)
    .fp <- function(p) round(p * 100, digits = 1)
    .fnp <- function(p) round(p, digits = 2)
    data.frame(PPV = .fp(tbl[2, 2]/sum(tbl[2, ])), PPV.CI = paste(sapply(binom.test(tbl[2, 
        2], sum(tbl[2, ]))$conf.int, .fp), collapse = "-"), NPV = .fp(tbl[1, 1]/sum(tbl[1, 
        ])), NPV.CI = paste(sapply(binom.test(tbl[1, 1], sum(tbl[1, ]))$conf.int, 
        .fp), collapse = "-"), FPR = .fp(tbl[2, 1]/sum(tbl[, 1])), FPR.CI = paste(sapply(binom.test(tbl[2, 
        1], sum(tbl[, 1]))$conf.int, .fp), collapse = "-"), FNR = .fp(tbl[1, 2]/sum(tbl[, 
        2])), FNR.CI = paste(sapply(binom.test(tbl[1, 2], sum(tbl[, 2]))$conf.int, 
        .fp), collapse = "-"), OR = .fnp(f$estimate), OR.CI = paste(sapply(f$conf.int, 
        .fnp), collapse = "-"), stringsAsFactors = FALSE)
    
}

.p2logitStage <- function(esets, fits, y = "debulking") {
    .doIt <- function(i) {
        gfits <- lapply(esets[-i], function(eset) {
            risk <- predict(fits[[i]]$fit, newdata = t(exprs(eset)))@lp
            glm(eset[[y]] ~ risk + tumorstage, eset, family = "binomial")
        })
        intercept.se <- do.call(rbind, lapply(gfits, function(x) summary(x)$coefficients[1, 
            ]))
        risk.se <- do.call(rbind, lapply(gfits, function(x) summary(x)$coefficients[2, 
            ]))
        idx <- sapply(gfits, function(x) nrow(summary(x)$coefficients) == 3)
        T.se <- do.call(rbind, lapply(gfits[idx], function(x) summary(x)$coefficients[3, 
            ]))
        res <- list(gfits, intercept = metafor::rma(yi = intercept.se[, 1], sei = intercept.se[, 
            2])$b, br = metafor::rma(yi = risk.se[, 1], sei = risk.se[, 2])$b, bt = metafor::rma(yi = T.se[, 
            1], sei = T.se[, 2])$b)
        
        #
        # predict(fit,newdata=data.frame(risk=fits[[i]]$risk@lp,stage=esets[[i]]$tumorstage),type='response')
        r <- fits[[i]]$risk@lp * res$br + esets[[i]]$tumorstage * res$bt + res$intercept
        exp(r)/(exp(r) + 1)
    }
    lapply(1:length(esets), .doIt)
}

.p2logitStage2 <- function(esets, fits, y = "debulking") {
    .doIt <- function(i) {
        esets.c <- .combineEsets(esets[-i], y = y)
        tmp.df <- data.frame(risk = predict(fits[[i]]$fit, newdata = esets.c$X)@lp, 
            covar = esets.c$y, tumorstage = as.factor(unlist(sapply(esets[-i], function(X) X$tumorstage))))
        fit <- glm(covar ~ risk + tumorstage, tmp.df, family = "binomial")
        predict(fit, newdata = data.frame(risk = fits[[i]]$risk@lp, tumorstage = as.factor(esets[[i]]$tumorstage)), 
            type = "response")
    }
    lapply(1:length(esets), .doIt)
}


.reclassPlotModel <- function(eset, m1, m2, ...) {
    r1 <- predict(m1, newdata = eset, type = "lp")@lp
    r2 <- predict(m2, newdata = eset, type = "lp")@lp
    .reclassPlot(eset$y, r1, r2, ...)
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
    panel.num = "A", ...) {
    labels <- sub(",.*20", " 20", sapply(dataset_ids, function(i) gsub("Cancer GenomeAtlas Research Network", 
        "TCGA", experimentData(esets[[i]])@lab)))
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
    
    par(mar = c(12.5, 5.5, 2, 0.5))
    plot(ids, unlist(x[1, ]), type = "l", las = 1, ylim = c(0.56, 0.65), xaxt = "n", 
        ylab = "", xlab = "", cex.axis = 1.2, cex.lab = 1.3, ...)
    axis(side = 1, at = ids, labels = labels[ids], las = 2, cex.axis = 1.2)
    title(ylab = ylab, cex.lab = 1.3, line = 4)
    points(ids, lb, col = "grey", lty = 2, type = "l")
    points(ids, ub, col = "grey", lty = 2, type = "l")
    mtext(panel.num, side = 3, cex = 1.8, at = 0)
    
}

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

.getRiskScore <- function(i, idx, esets, fits) {
    # X = t(exprs(esets[[i]]))
    predict(fits[[idx]]$fit, newdata = esets[[i]], type = "lp")@lp
}

addClinical <- function(idx, esets, fits) {
    df <- .createClinical(esets[-idx])
    df$risk <- unlist(sapply((1:length(esets))[-idx], .getRiskScore, idx, esets, 
        fits))
    fit <- coxph(y ~ tumorstage + debulking + 
        risk, data = df)
    dfnew <- .createClinical(esets[idx])
    dfnew$risk <- fits[[idx]]$risk@lp
    list(fit = fit, risk = predict(fit, newdata = dfnew))
}

addClinicalBinary <- function(idx, esets, fits) {
    df <- .createClinical(esets[-idx], clinical.covars = c("age_at_initial_pathologic_diagnosis", 
        "debulking", "tumorstage", "days_to_death", "vital_status", "oneyr"))
    df$risk <- unlist(sapply((1:length(esets))[-idx], .getRiskScore, idx, esets, 
        fits))
    fit <- glm(as.logical(oneyr) ~ as.factor(age_at_initial_pathologic_diagnosis > 
        60) + debulking + risk, data = df, family = "binomial")
}

addClinicalBinaryBK <- function(idx, esets, fits) {
    df <- .createClinical(esets[-idx], clinical.covars = c("age_at_initial_pathologic_diagnosis", 
        "debulking", "tumorstage", "days_to_death", "vital_status", "oneyr"))
    df$risk <- unlist(sapply((1:length(esets))[-idx], .getRiskScore, idx, esets, 
        fits))
    fit <- glm(df$debulking == "suboptimal" ~ as.factor(age_at_initial_pathologic_diagnosis > 
        60) + debulking + risk, data = df, family = "binomial")
}

addClinicalAll <- function(esets, fits) {
    res <- lapply(1:length(esets), addClinical, esets, fits)
    esets.r <- lapply(1:length(esets), function(i) {
        esets[[i]]$risk <- res[[i]]$risk
        esets[[i]]
    })
    names(esets.r) <- names(esets)
    esets.r
}

addClinicalAllBinary <- function(esets, fits) {
    res <- lapply(1:length(esets), addClinicalBinary, esets, fits)
    esets.r <- lapply(1:length(esets), function(i) {
        esets[[i]]$risk <- res[[i]]$risk
        esets[[i]]
    })
    names(esets.r) <- names(esets)
    esets.r
}

.pOverlap <- function(A, B, background) {
    n1 <- length(A)
    n2 <- length(B)
    n <- length(background)
    m <- length(intersect(A, B))
    phyper(min(n1, n2), n1, n - n1, n2) - phyper(m - 1, n1, n - n1, n2)
}
 
