## explore downloaded storet data to check on
## cond measurements

explorestoret <- function() {
    print(table(instoret$CharacteristicName))

    incvec <- instoret$CharacteristicName == "Specific conductance"
    df1 <- instoret[incvec,]

    print(table(df1$ResultMeasure.MeasureUnitCode))
    incvec <- df1$ResultMeasure.MeasureUnitCode == "uS/cm"
    print(sum(incvec))
    df1 <- df1[incvec,]

    print(summary(as.numeric(df1$ResultMeasureValue)))

    df1$val <- as.numeric(df1$ResultMeasureValue)
    incvec <- df1$val < 10
    incvec[is.na(incvec)] <- F
    df1$val[incvec] <- NA

    hist(log(df1$val), breaks = 40)
    print(median(df1$val, na.rm = T))
}

explorestoret()
