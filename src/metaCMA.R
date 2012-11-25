library(metafor)
library(survHD)
library(survcomp)
library(rmeta)

.combineEsets <- function(esets, y="y",
    probesets=featureNames(esets[[1]]),ComBat=TRUE) {
    if (length(esets)==0) return(NULL)
    X = do.call("cbind", lapply(esets, exprs))
    if (!is.null(y) && class(esets[[1]][[y]])=="Surv") {
        y <- do.call("rbind", lapply(esets, function(x) x[[y]]))
        y <- Surv(y[,1], y[,2])
    } else {
        y <- do.call("c", lapply(esets, function(x) as.character(x[[y]])))
    }
    batch =  as.factor(do.call(c,lapply(1:length(esets), function(i)
        rep(names(esets)[i], ncol(esets[[i]])))))

    l = list(X=t(X[probesets,]),y=y,batch=batch)
    colnames(l$X) = make.names(colnames(l$X), unique = FALSE)
    l
}

.calcRMA <- function(idx, coefs, rma.method="FE") {
    lapply(1:nrow(coefs$c),function(i) try(metafor::rma(yi = coefs$c[i,idx], sei =
        coefs$se[i,idx], method=rma.method)))
}

.createXy <-function(idx, esets, y="y", coefs, n,
rma.method="FE", filter.fun=.defaultFilter, modeltype="compoundcovariate", only.sign=NULL){ 
    probesets = featureNames(esets[[1]])
    coefficients = NULL
    model = NULL
    pvalues <- NULL

    idx.not.want = -which(sapply(esets,filter.fun))
    # idx is minus validation set
    idx.t = c(idx, idx.not.want)
    if (length(esets[idx.t]) == 0) {
        warning("No training data passed the filter.")
        return(NULL)
    }
    idx.v = c()
    if (length(idx) > 0) idx.v = -idx
    
    if (!is.null(n)) {
        res.rma = .calcRMA(idx.t, coefs, rma.method)
        pvalues = as.numeric(sapply(res.rma, function(r) try(r$pval)))
        #ids = apply(coefs$c[,idx],1,function(x) sum(sign(x))==length(x) || sum(sign(x))==-length(x) )
        names(pvalues) = rownames(coefs$c)
        coefficients = sapply(res.rma, function(x) x$b)
        #pvalues[!ids] = 1
        names(coefficients) = make.names(rownames(coefs$c))
        if (!is.null(only.sign)) {
            cat("Usign only", only.sign, "signed coefficients.")
            idx = sign(coefficients) == sign(only.sign)
            pvalues = pvalues[idx]
            coefficients = coefficients[idx]
        }
        probesets = head(order(pvalues),n)
        model = new("linearriskscore",
        coefficients=coefficients[probesets],modeltype=modeltype)
    } else {
        cat("No RMA because no n")
    }
    xy = list(train    = .combineEsets(esets[idx.t], y, probesets),
              validate = .combineEsets(esets[idx.v], y, probesets),
              model    = model, idx.t = idx.t, idx.v = idx.v, 
              pvalues = pvalues
    )
    xy
}

metaCMA.train <- function(i, esets, y="y", coefs, n, method=NULL, rma.method="FE",
verbose=FALSE, filter.fun=.defaultFilter, modeltype="compoundcovariate", only.sign=NULL,...) {
    if (verbose) cat("Dataset", i, "of", length(esets),"...\n")
    tmp <- .createXy(-i,  esets, y, coefs, n, rma.method,
    filter.fun,modeltype=modeltype, only.sign=only.sign)
    if (is.null(tmp)) {
        return(list)
    } else if (is.null(method)) {
        fit = tmp$model
    } else {
        fit = do.call(method,args=list(X=tmp$train$X, y=tmp$train$y, ... ))
    }
    list(fit = fit, risk = predict(fit, tmp$validate$X, type="lp"))
}

.defaultCoef <- function(x, y, X,...) {
    if (is.Surv(y)) { 
        fit = coxph(y~x) 
        return(summary(fit)$coefficients[c(1,3)])
    }    
    fit = glm(y~x,...)
    cf = summary(fit)$coefficients
    if (nrow(cf)==2) return(cf[2,1:2])
    warning(paste("Coefficients not found",cf))
    c(NA,NA)
}

metaCMA.coefs <- function(esets, y="y", coef.fun=.defaultCoef, ...) {
    res = lapply(esets, function(X) apply(exprs(X), 1, coef.fun, X[[y]], X, ...))
    coefs = list(c = mapply(cbind, lapply(res, function(x) x[1,])),
                 se = mapply(cbind, lapply(res, function(x) x[2,])))
    rownames(coefs$c)  = featureNames(esets[[1]])
    rownames(coefs$se) = featureNames(esets[[1]])
    # only probe sets without missing values and non-zero se in all datasets
    idx = complete.cases(coefs$c) & apply(coefs$se,1,function(x) sum(x==0)==0)
    coefs$c  = coefs$c[idx,]
    coefs$se = coefs$se[idx,]
    coefs
}

.defaultFilter <- function(eset) {
    ncol(eset) < 75
}

metaCMA.opt <- function(esets, ...) {
    .createXy(esets,idx=c(),...)
}

metaCMA <- function(esets, y="y",  coefs=NULL, n, method=NULL,
rma.method="FE", filter.fun=.defaultFilter, modeltype="compoundcovariate",
only.sign=NULL, ... ) {
    if (is.null(coefs)) coefs = metaCMA.coefs(esets, y)
    fits  = lapply(1:length(esets), metaCMA.train, esets=esets, y=y, coefs=coefs,
    n=n, method=method,
    rma.method=rma.method,filter.fun=filter.fun, modeltype=modeltype,
    only.sign=only.sign, ...)
    names(fits) = names(esets)
    list(fits=fits, y=y, rma.method=rma.method)
}

metaCMA.eval <- function(ids, esets, y="y", object, tau=365.25*4,...) {
    lapply(ids, function(i)
    evaluate(object$fits[[i]]$risk, measure=new("UnoC"),
        newy=esets[[i]][[y]], add=list(tau=tau) ))
}

.camelCase <- function(x) paste(toupper(substr(x,1,1)), substr(x,2,nchar(x)), sep="")

.forestplot <- function(esets, y, label, rma.method="FE", measure="hr",
at = NULL, concordance=TRUE, xlab=ifelse(concordance, "Concordance", "Hazard ratio"),...) {
    if (is.null(names(esets))) names(esets) = paste("Study", 1:length(esets))
    names(esets) = gsub("_", " ", names(esets))
    if (concordance) {
        coefs = sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~esets[[i]][[label]]))$concordance)   
        res.rma = metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        if (is.null(at)) at = seq(0.4,1.0,0.1)
        forest.rma(res.rma, 
        xlab=xlab,  slab=sapply(names(esets), .camelCase),
        refline=0.5, at=at,...)
        return(res.rma)
    } else {
        coefs = sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~esets[[i]][[label]]))$coefficients[c(1,3)])   
        res.rma = metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        if (is.null(at)) at = log(c(0.25,1,4,20))
        forest.rma(res.rma, 
        xlab=xlab,  slab=sapply(names(esets), .camelCase),
        atransf=exp, at=at, ...)
        return(res.rma)
    }
}

metaCMA.concordance <- function(esets, y="y", risks, rma.method="FE") {
        coefs <- sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~risks[[i]]))$concordance)
        res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        list(res.rma, coefs)
}

metaCMA.hr <- function(esets, y="y", risks, rma.method="FE") {
        coefs <- sapply(1:length(esets), function(i)
        summary(coxph(esets[[i]][[y]]~risks[[i]]))$coefficients[c(1,3)])
        res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,], method=rma.method)
        list(res.rma, coefs)
}

metaCMA.forest <- function(esets, metacma, y="y", mlab="Overall", ...) {
    tmp = names(esets)
    esets = lapply(1:length(esets), function(i) { esets[[i]]$risk =
        metacma$fits[[i]]$risk@lp; esets[[i]]} )   
    names(esets) = tmp    
    .forestplot(esets, y, label="risk", rma.method=metacma$rma.method, mlab=mlab, ...)
}

metaCMA.forest.models <- function(esets, y="y", risks1, risks2,
mlab="Overall",concordance=TRUE,...) {
    tmp <- names(esets)
    labeltext <- cbind(c("", sapply(tmp, function(x)
    c(x,NA,NA)),"Overall",NA))
                     #  c("Signature", rep(c("Meta-Analysis","TCGA",NA),
                     #  length(tmp)), c("Meta-Analysis", "TCGA") ))
    if (concordance) {
        rma1 <- metaCMA.concordance(esets,y, risks1)
        rma2 <- metaCMA.concordance(esets,y, risks2)
    } else {
        rma1 <- metaCMA.hr(esets,y, risks1)
        rma2 <- metaCMA.hr(esets,y, risks2)
    }
    r <- do.call(rbind, lapply(1:ncol(rma1[[2]]),function(i)
    rbind(rma1[[2]][,i],rma2[[2]][,i], c(NA,NA))))
    r.mean <- c(NA,r[,1], rma1[[1]]$b,  rma2[[1]]$b)
    r.lower <- c(NA,r[,1]-(r[,2]*1.96), rma1[[1]]$ci.lb,  rma2[[1]]$ci.lb)
    r.upper <- c(NA,r[,1]+(r[,2]*1.96), rma1[[1]]$ci.ub,  rma2[[1]]$ci.ub)
    col=meta.colors(line=c(rep(c(NA, "darkblue", "seagreen"),length(tmp)+1)), zero="firebrick", box=c(rep(c(NA," royalblue", "forestgreen"),length(tmp)+1)))

    forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower,
    upper=r.upper,zero=ifelse(concordance,0.5,0), col=col,
    xlog=concordance==FALSE,
    align=c("l"),  xlab=ifelse(concordance, "Concordance Index", "Hazard Ratio"),
    is.summary=(c(rep(FALSE,length(tmp)*3+1), TRUE, TRUE)),...)
}

metaCMA.forest.probeset <- function(esets, y="y",
probeset="208079_s_at", mlab="Overall", rma.method="FE",...) {
    esets = esets[sapply(esets, function(x) probeset %in% featureNames(x))]
    esets = lapply(esets, function(x) { x[[probeset]] = exprs(x)[probeset,]; x })

    .forestplot(esets, y, label=probeset, rma.method=rma.method, mlab=mlab, ...)
}

metaCMA.forest.label <- function(esets, y="y",
label, mlab="Overall", rma.method="FE",...) {
    .forestplot(esets, y, label, rma.method=rma.method, mlab=mlab, ...)
}

metaCMA.common.gene.esets <- function(esets) {
    fns = featureNames(esets[[1]])
    tmp = mapply(function(x) fns <<- intersect(fns,x), lapply(esets, featureNames))
    lapply(esets, function(x) x[fns,])
}

metaCMA.xtable <- function(xtbl, hline.after=getOption("xtable.hline.after", c(-1,0,nrow(xtbl))),
...) {
    print(xtbl,
                  hline.after=NULL,
                  add.to.row=list(pos=lapply(hline.after, function(x) x),
                  command=c("\\toprule\n",
                            rep("\\midrule\n", length(hline.after)-2),
                            "\\bottomrule\n")),...)
}

metaCMA.censor <- function(esets, y="y", censor.at=NULL) {
    .censor <- function(eset) {
        if (!is.null(censor.at)) { 
            eset[[y]][eset[[y]][,1] > censor.at,2] <- 0
            eset[[y]][eset[[y]][,1] > censor.at,1] <- censor.at
        }
        eset
    }
    lapply(esets,.censor)
}

# simple differential expression over all datasets using the limma package
metaCMA.limma <- function(esets, groups, contrasts) {
    require(limma)
    .doLimma <- function(eset) {
        TS <- as.factor(eset[[groups]])
        design <- model.matrix(~0+TS)
        colnames(design) <- levels(TS)
        fit <- lmFit(exprs(eset), design = design)
        cont <-
        makeContrasts(contrasts=contrasts,levels=design)
        fit2 <- contrasts.fit(fit, cont)
        fit2 <- eBayes(fit2)
    }    
    lapply(esets, .doLimma)
}

.getCoefsSubset <- function(coefs, idx) {
    ret <- coefs
    ret[[1]] <- ret[[1]][,idx]
    ret[[2]] <- ret[[2]][,idx]
    ret
}

metaCMA.powerset <- function(n) {
    unlist(lapply(lapply(1:n,function(i) combn(x=1:n, i)), function(x)
        lapply(1:ncol(x), function(i)  x[,i])),recursive=FALSE)
}

metaCMA.allcombinations <- function(i, esets, eval.fun = metaCMA.eval, ...) {
    .doPS <- function(ps) {
        idx <- c(i,(1:length(esets))[-i][ps])
        tmp <- metaCMA.train(1, esets[idx], coefs=.getCoefsSubset(coefs,idx),filter.fun=function(eset)
        return(FALSE),... )
        eval.fun(1,esets[idx],object=list(fits=list(tmp)))
     }
     pss <- metaCMA.powerset(length(esets[-i]))
     ret <- lapply(pss, .doPS)
     n <- lapply(pss, function(ps) sum(sapply(esets[-i][ps], ncol)))
     list(evaluation=ret, n=n)
}

