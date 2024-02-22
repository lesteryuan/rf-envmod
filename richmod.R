## explore biological metrics values
## questions for IN:
## How is biological data used currently in water programs?
exploremet <- function() {
    ## figure out unique sample id
    richness.in$idnew <- paste(richness.in$lsite, richness.in$date, sep = "-")

    ## need to ask mitchell about repeat samples on the same day

    ## pick first sample
    richness.in <- richness.in[!duplicated(richness.in$idnew),]

    metrics$idnew <- paste(metrics$lsite, metrics$date, sep = "-")
    nsamp <- table(metrics$idnew)
    idp <- names(nsamp)[nsamp > 1]
    incvec <- rep(F, times = nrow(metrics))

    for (i in idp) incvec <- incvec|metrics$idnew == i

    df1 <- metrics[incvec,]
    print(summary(tapply(df1$individuals.count, df1$idnew, function(x) diff(range(x)))))
    stop()
    incvec <- metrics$idnew == idp[1]
    print(metrics[incvec,])

    plot(log(metrics$individuals.count), log(metrics$taxa.count))
    stop()

}
#exploremet()

## model relationship between richness and environmental
## variables

richmod <- function(envdat.tss, envdat.cond, pred, richdat) {

    varpred <- names(pred)[-c(1,2,4)]

    ## merge data
    trich <- tapply(richdat$TAXA_COUNT, richdat$LSite, mean)
    pintol <- tapply(metrics$pct.intolerant, metrics$lsite, mean)

    dftemp <- data.frame(LSite = names(trich), trich=as.vector(trich))
    print(nrow(dftemp))
    dftemp <- merge(dftemp, data.frame(LSite = names(pintol),
                                       pintol = as.vector(pintol)),
                    by = "LSite")
    
    df1 <- merge(envdat.tss[, c("LSite", "tss.pred")],dftemp, 
                 by = "LSite")
    print(nrow(df1))
    df1 <- merge(df1, envdat.cond[, c("LSite", "cond.pred")], by = "LSite")
    print(nrow(df1))
    df1 <- merge(df1, unique.data.frame(pred[, varpred]), by = "LSite")
    print(nrow(df1))

    df2 <- df1[, -1]
    df2 <- df2[, -which("kffactws" == names(df2))]
    df2 <- df2[, -which("pctsilicicws" == names(df2))]
    print(summary(df2))

    pvarnames <- names(df2)
    pvarnames <- pvarnames[-which(pvarnames == "pintol")]
    pvarnames <- pvarnames[-which(pvarnames == "trich")]
    print(pvarnames)

    incvec <- df2$LATITUDE_MEASURE < 39.6
    dflist <- split(df2, incvec)
    ## North is [[1]], South is [[2]]
    subsetname <- c("North", "South")

    npred <- c(6,10)
    ## names of anthropogenic variables to plot for each subset
#    varplot <- list(c("inorgnwetdep.2008ws", "pctcrop2019ws",
#                      "pctnonagintrodmanagvegws"),
#                    c("tss.pred", "pctconif2019ws"))
#    lab0 <- list(c("Wet N deposition", "Percent crop",
#                   "Percent introduced managed vegetation"),
#                 c("Predicted TSS", "Percent coniferous forest"))
    varplot <- list(c("pctcrop2019ws","cond.pred", "tss.pred"),
                    c("tss.pred", "cond.pred"))
    lab0 <- list(c("Percent crop", "Conductivity",
                   "TSS"),
                 c("TSS", "Conductivity"))
    require(pdp)
    dev.new()
    plot(dflist[[1]]$cond.pred, dflist[[1]]$pintol, axes = F)
    logtick.exp(0.001, 10, c(1), c(T,F))
    abline(h = c(16, 32))
    axis(2)
    stop()

    dev.new()
#    par(mar = c(4,4,2,1), mfrow = c(1,2), mgp = c(2.3,1,0),
#        bty = "l")
    par(mar = c(4,4,1,1), mfrow = c(1,3), mgp = c(2.3,1,0),
        bty = "l")
    for (i in 1:1) {
#    for (i in 1:length(dflist)) {
        set.seed(10)
        require(ranger)
        mod <- ranger(data = dflist[[i]][, c(pvarnames, "pintol")],
                      dependent.variable.name = "pintol",
                      num.trees = 5000, importance = "permutation")
        print(mod)

        plot(mod$predictions, dflist[[1]]$pintol)
        abline(0,1)

        err0 <- mod$prediction.error
        
        varsel <- names(rev(sort(mod$variable.importance)))
        print(varsel)
        stop()
        print(varsel[1:npred[i]])

        ## set to true to find the number of predictors
        ## that achieves the same prediction error as the full model
        findnvar <- FALSE
        if (findnvar) {
            nvar <- 5:16
            err.sc <- rep(NA, times = length(nvar))
            for (j in 1:length(nvar)) {
                mod <- ranger(data = dflist[[i]][, c(varsel[1:nvar[j]],
                                  "pintol")],
                              dependent.variable.name = "pintol",
                              num.trees = 5000, importance = "permutation")
                err.sc[j] <- mod$prediction.error/err0
            }
            plot(nvar, err.sc)
            stop()
        }
        else {
            ## number of predictors is known and saved in npred
            ## so refit the model with npred
            mod <- ranger(data = dflist[[i]][, c(varsel[1:npred[i]],
                                  "pintol")],
                          dependent.variable.name = "pintol",
                          num.trees = 5000, importance = "permutation")
            print(mod)
        }
        ## set plotpred to TRUE to plot predicted vs. observed
        ## otherwise, will plot partial dependence relationships
        ## for anthropogenic variables
        plotpred <- F
        if (plotpred) {
            plot(mod$predictions, dflist[[i]]$pintol,
                 xlab = "Richness (predicted)", ylab = "Richness (observed)",
                 pch = 21, col = "grey", bg = "white")
            abline(0,1, lty = "dashed")
            mtext(subsetname[i], side = 3, line = 0)

        }
        else {
            for (k in 1:length(varplot[[i]])) {
                pout <- partial(mod, pred.var = varplot[[i]][k], plot = F,
                                approx = T)
                print(names(pout))
                plot(pout[,1], pout[,2] - mean(pout[,2]), type = "l",
                     xlab = lab0[[i]][k], ylim = c(-5, 7),ylab = "Change in percent intolerant")
            }
        }
    }

    stop()
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(mod$predictions, df2$TAXA_COUNT, col = "grey")
    abline(0,1)


    stop()

    plot(df1$cond.obs, df1$tss.obs)
    mod <- gam(TAXA_COUNT ~ s(tss.pred, k = 4) + s(cond.pred, k = 4),
               data = df1, family = "poisson")
               print(summary(mod))

    stop()

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0), bty = "l")
    plot(df1$tss.obs, df1$TAXA_COUNT, col = "grey")

    require(mgcv)
    mod <- gam(TAXA_COUNT ~ s(tss.obs, k = 4), data = df1,
               family = "poisson")
    xnew <- seq(min(df1$tss.obs, na.rm = T), max(df1$tss.obs, na.rm = T),
                length = 50)

    predout <- predict(mod, list(tss.obs = xnew), type = "response")
    lines(xnew, predout)
    print(summary(mod))

    plot(df1$tss.pred, df1$TAXA_COUNT, col = "grey")

}

richmod(envdat.tss = dftss, envdat.cond = dfcond, pred = pred,
        richdat = richness.in)
