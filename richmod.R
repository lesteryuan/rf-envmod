## model relationship between richness and environmental
## variables

richmod <- function(envdat.tss, envdat.cond, pred, richdat) {

    varpred <- names(pred)[-c(1,2,4)]

    ## merge data
    df1 <- merge(envdat.tss[, c("LSite", "tss.pred")],
                 richdat[, c("LSite", "TAXA_COUNT")], by = "LSite")
    df1 <- merge(df1, envdat.cond[, c("LSite", "cond.pred")], by = "LSite")
    df1 <- merge(df1, pred[, varpred], by = "LSite")

    df2 <- df1[, -1]
    df2 <- df2[, -which("kffactws" == names(df2))]
    df2 <- df2[, -which("pctsilicicws" == names(df2))]
    print(summary(df2))

    incvec <- df2$LATITUDE_MEASURE < 39.6
    dflist <- split(df2, incvec)
    ## North is [[1]], South is [[2]]
    subsetname <- c("North", "South")

    npred <- c(10,10)
    ## names of anthropogenic variables to plot for each subset
    varplot <- list(c("inorgnwetdep.2008ws", "pctcrop2019ws",
                      "pctnonagintrodmanagvegws"),
                    c("tss.pred", "pctconif2019ws"))
    lab0 <- list(c("Wet N deposition", "Percent crop",
                   "Percent introduced managed vegetation"),
                 c("Predicted TSS", "Percent coniferous forest"))
    require(pdp)

    dev.new()
#    par(mar = c(4,4,2,1), mfrow = c(1,2), mgp = c(2.3,1,0),
#        bty = "l")
    par(mar = c(4,4,1,1), mfrow = c(1,3), mgp = c(2.3,1,0),
        bty = "l")
    for (i in 2:2) {
#    for (i in 1:length(dflist)) {
        set.seed(10)
        require(ranger)
        mod <- ranger(data = dflist[[i]],
                      dependent.variable.name = "TAXA_COUNT",
                      num.trees = 5000, importance = "permutation")
        print(mod)
        err0 <- mod$prediction.error
        varsel <- names(rev(sort(mod$variable.importance)))

        print(varsel[1:npred[i]])
        ## set to true to find the number of predictors
        ## that achieves the same prediction error as the full model
        findnvar <- FALSE
        if (findnvar) {
            nvar <- 5:12
            err.sc <- rep(NA, times = length(nvar))
            for (j in 1:length(nvar)) {
                mod <- ranger(data = dflist[[i]][, c(varsel[1:nvar[j]],
                                  "TAXA_COUNT")],
                              dependent.variable.name = "TAXA_COUNT",
                              num.trees = 5000, importance = "permutation")
                err.sc[j] <- mod$prediction.error/err0
            }
            plot(nvar, err.sc)
        }
        else {
            ## number of predictors is known and saved in npred
            ## so refit the model with npred
            mod <- ranger(data = dflist[[i]][, c(varsel[1:npred[i]],
                                  "TAXA_COUNT")],
                          dependent.variable.name = "TAXA_COUNT",
                          num.trees = 5000, importance = "permutation")
            print(mod)
        }
        ## set plotpred to TRUE to plot predicted vs. observed
        ## otherwise, will plot partial dependence relationships
        ## for anthropogenic variables
        plotpred <- F
        if (plotpred) {
            plot(mod$predictions, dflist[[i]]$TAXA_COUNT,
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
                plot(pout[,1], pout[,2], type = "l",
                     xlab = lab0[[i]][k], ylab = "Change in richness")
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
