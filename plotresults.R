

## test site results 12.19.2023
plotresults2 <- function() {

    obsmn <- tapply(condpred$Conductance..uS.cm.,
                    condpred$STATION_NAME, function(x) mean(log(x)))
    predmn <- tapply(condpred$cond_test,
                     condpred$STATION_NAME, mean)


    dftemp <- data.frame(LSite = names(obsmn),
                         obsmn = as.vector(obsmn), predmn = as.vector(predmn))
    dftemp$delt <- dftemp$obsmn - dftemp$predmn
    print(nrow(dftemp))

    latmn <- tapply(locs$LATITUDE_MEASURE, locs$LSite, mean)
    lonmn <- tapply(locs$LONGITUDE_MEASURE, locs$LSite, mean)
    locsmn <- data.frame(LSite = names(latmn), latmn = as.vector(latmn),
                         lonmn = as.vector(lonmn))
    dftemp <- merge(dftemp, locsmn, by = "LSite")
    print(nrow(dftemp))
    print(summary(dftemp$delt))
    incvec <- dftemp$delt > 0.75
    incvec2 <- dftemp$delt < -0.75

    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    plot(dftemp$predmn, dftemp$obsmn, pch = 21, col = "grey39", bg = "white", axes = F,
         xlab = "Predicted conductivity", ylab = "Observed conductivity")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1, lty = "dashed")
    points(dftemp$predmn[incvec], dftemp$obsmn[incvec], pch = 16, col = "red")
    points(dftemp$predmn[incvec2], dftemp$obsmn[incvec2], pch = 16, col = "green4")

    require(maps)
    require(mapproj)
    dev.new()
    par(mar = c(1,1,1,1))
    map("state", region = "indiana", proj = "albers", par = c(30,40))
    pout <- mapproject(dftemp$lonmn, dftemp$latmn, proj = "")
    points(pout$x, pout$y, pch = 21, col = "grey", bg = "white")
    points(pout$x[incvec], pout$y[incvec], pch =16, col = "red")
    points(pout$x[incvec2], pout$y[incvec2], pch =16, col = "green4")
    return()
}
plotresults2()

## plot environmental model results
plotresults <- function() {

#    mcodat <- read.csv("MCO 2022 sites model test.csv")
    mcodat <- read.csv("MCO ref22 model test.csv")
    print(names(mcodat))

    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot( log(envdata1819$tss.result), log(envdata1819$ptl),
         pch = 21, col = "grey", bg = "white", axes = F,
         xlab = "TSS (mg/L)", ylab = expression(Total~P~(mu*g/L)))
    points( log(allchem$tss),log(allchem$tp*1000), col = "red")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    stop()

#    png(width = 6, height = 3, pointsize = 10, units = "in",
#        res = 600, file = "pred.obs.png")
    dev.new()
    par(mar = c(4,4,1,1), mfrow= c(1,2), mgp = c(2.3,1,0))
    plot(mcodat$tss.result_test, log(mcodat$Avg.TSS), axes = F,
         pch = 21, col = "grey39", bg = "white",
         xlab ="TSS (predicted)",ylab = "TSS (observed)")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1)
    plot(mcodat$cond_test, log(mcodat$Avg.SC), axes = F,
         pch = 21, col = "grey39", bg = "white",
         xlab = "Cond (predicted)", ylab = "Cond (observed)")
    logtick.exp(0.001, 10, c(1,2), c(T,T))
    abline(0,1)
    stop()
    dev.off()

}

#plotresults()
