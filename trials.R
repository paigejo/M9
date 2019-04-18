# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = TRUE, 
#                       G=G, debugPlotting=TRUE, logPenaltyPar=log(1), logDiffPenaltyPar=log(1), 
#                       sharedSpatialProcess=TRUE, jointShared = TRUE)
# this one works well in that sd is small relative to mean, but has biased GPS predictions (all dStar=25000:
# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE, 
#                       G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
#                       sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE, 
#                       doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar, 
#                       dStarGPS=dStarGPS)
# nKnotsMeanGPS=7 accidentally: sd wayyyy too small, GPS mean very small, x residuals scrunch around 0 but aren't biased, y residuals unbiased but heteroscedastic:
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, nKnotsMeanGPS=7)
# same as before but nKnotsMeanGPS=1. 
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE)
# same as before but with a variance spline: still reasonable
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE)
# same as before but including locking rate variance inflation
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, inflateVarLocking=TRUE)
# same as before but fixed taper penalty and including different variance for gps data: now gps data is almost always predicted as very close to 0...
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE)
# allow for mean spline: biased x residuals, with bias increasing towards the north, and x residuals front greatly toward 0
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE)
# include gamma spline, increase number of knots: there is bias (decreasing x residuals as fun of lat) causing issues, poor x predictions
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=11, nKnotsVarGPS=7, nKnotsGamma=7)
# get rid of the different GPS variance spline: 
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=FALSE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=11, nKnotsGamma=7)
# increase number of GPS knots: looks like smoothness penalty should be increased.  X residuals downardly biased, shrunk to 0
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=FALSE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=15, nKnotsGamma=15)
# try only including gamma spline, different variance, increase smoothness penalty slightly: same as last time
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=7, nKnotsGamma=15)
# try not using a penalty: results turn non-physical, y residuals look good, but x residuals are shifted up and somewhat shrunk towards 0
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=7, nKnotsGamma=15, doSmoothnessPenalty=FALSE)
# now try an unpenalized simple model: mean nonphysical in the south, y residuals look good, x residuals still have structure and shrunk somewhat toward 0
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=FALSE)
# now allowing for shared spatial structure: it exploded...
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = FALSE, includeGammaSpline = FALSE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=FALSE)
# now estimate the proportion of shared variance in the gps data, remove mean spline (keep variance spline): explosion...
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=FALSE)
# add in smoothness penalty: still exploded, all splines ~constant
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=TRUE)
# add mean spline, keep smoothness penalty
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = FALSE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=TRUE)
# add in smoothness penalty: still exploded, all constant
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = FALSE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=TRUE)
# add in mean spline: unidentifiable?? outer mgc went to infinity...
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = FALSE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=1, doSmoothnessPenalty=TRUE)
# remove mean spline, add in conditional estimation for gamma spline: Results look great for the most part, although x residuals slightly shrunk toward 0 and locking rate inflation is 15.84. mean and standard deviation appear constant, proportion of shared variance is roughly 1
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# stop estimating locking variance inflation, remove smoothness penalty and use same taper, add in mean spline: nonphysical result
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=FALSE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=FALSE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=FALSE, conditionalGamma=TRUE)
# add slight smoothness penalty: still terrible
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=FALSE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=FALSE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# go back a couple of steps to when it was working, but stop estimating locking variance inflation: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=FALSE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# try adding a difference penalty: didn't help
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = TRUE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.01),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=FALSE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# go back to the working model but make the variance shared: residuals are heteroscedastic and biased as a function of lat. Shrunk slightly towards 0, and discovered inflateVarLocking is the variance inflation, not SD inflation so fixed code. var inflation was almost 9. SD larger than mean in south
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=FALSE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# try previous model that worked (same as above but with different variance splines), and increase number of variance splines: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.01), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# increase smoothness penalty: splines are now flat, but results look good otherwise. residuals heteroscedastic and x residuals shrunk somewhat towards 0
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE)
# decrease smoothness penalty just a bit, try not constraining the mean: mean goes to 0, sd looks reasonable, very biased y residuals
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE)
# constrain the mean: results again look pretty good, although x residuals somewhat shrunk to 0, and some heteroscedasticity
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# decreased smoothness penalty slightly
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.05), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# decreased smoothness penalty slightly again
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# decreased smoothness penalty slightly again: pretty much the same
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# remove smoothness penalty: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=FALSE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# add smoothness penalty back in but smaller: looks basically the same still
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.00001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# decrease smoothness penalty
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.00000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
testFit = fullFit
fullFit$obj$theta = fullFit$opt$par
fullFit$obj$par = fullFit$opt$par
testFit = fitModelTMB(fullFit, fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.00000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=8)
# try increasing the smoothness penalty and the number of variance GPS knots: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# increase smoothness penalty more
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(100), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE)
# allow only 1 taper parameter for each dataset, but allow variance, gamma splines for both: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(100), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=TRUE, estimateGpsShared=TRUE, includeGammaSpline=TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, nKnotsGPS=1, nKnots=1)
# don't use seperate taper for GPS data: nonphysical
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(100), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=TRUE, estimateGpsShared=TRUE, includeGammaSpline=TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, nKnotsGPS=1, nKnots=1)
# allow a seperate taper but only shifted form the original: whoops, too few variance knots for the fault.  Looks alright otherwise
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=TRUE, estimateGpsShared=TRUE, includeGammaSpline=TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, nKnotsGPS=1, nKnots=5)
# fix the variance knot issue above: results look solid except everything but the taper is constant anyway.  Also the fault sd and mean are equal
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=TRUE, estimateGpsShared=TRUE, includeGammaSpline=TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsVar=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, nKnotsGPS=1, nKnots=5)
# try to see what's going on with the penalty for the best model: looks fine, penalty: 0.03641983
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.00000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, maxCount=1)
# reduce the penalty: penalty term: 0.02550234
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0000000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, maxCount=1)
# reduce the penalty more: penalty term now 0.01261327
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.000000000000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, maxCount=1)
# remove the variance spline smoothness penalty: penalty term now 0.01878
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.000000000000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=TRUE, maxCount=1, varSmoothnessPenalty=FALSE)
# try removing constraint: it works!!!
fullFit$obj$par = fullFit$opt$par
fullFit$obj$theta = fullFit$opt$par
fullFit = fitModelTMB(fullFit, fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.000000000000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE)
# try removing variance knots and smoothness penalty: still looks great
initialParameters = fullFit$opt$par
initialParameters = initialParameters[-(13:16)]
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.000000000000001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=FALSE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=1, paramInit=initialParameters)
# increase smoothness penalty: nonphysical
initialParameters = fullFit$opt$par
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(5), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=1, paramInit=initialParameters, 
                      recompile=TRUE)
# decrease smoothness penalty: now SDs are ~30 (but it doesn't matter since shared term included).  Also, SEs are defined except for logitOmega!
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=1, paramInit=initialParameters, 
                      recompile=TRUE)
# Allow for variance spline
initialParameters = fullFit$opt$par
# initialParameters = c(initialParameters[1:12], "betasd"=rep(0, 4), initialParameters[13:20])
initialParameters[13:16] = c(-2.374877e-11, -1.341930e-10, -6.077585e-11, 2.359263e-10)
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=FALSE)
# allow variance spline for gps data
initialParameters = fullFit$opt$par
initialParameters = c(initialParameters[1:17], "betasdGPS"=rep(0, 4), initialParameters[18:24])
names(initialParameters)[18:21] = "betasdGPS"
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=5, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=TRUE)
# make variance be constant for gps, but flexible for the subsidence (reparameterize the variance), remove smoothness penalty for variance: variance is still constant...
initialParameters[17] = initialParameters[17] + initialParameters[12]
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=TRUE, reparameterizeVar=TRUE, varLogLambda=log(.0001))
# restart the optimization based on these initial parameters modified slightly: very large standard deviations, but physical mean
initialParameters = fullFit$opt$par - c(rep(0, 12), 0.2*fullFit$grad[13:16], rep(0, 8))
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=FALSE, reparameterizeVar=TRUE, varLogLambda=log(.0001))
# add-in slight variance smoothness penalty: standard deviation now varies slightly!
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=TRUE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=TRUE, reparameterizeVar=TRUE, varLogLambda=log(.0001))
# try reducing the variance smoothness penalty: results look basically the same
initialParameters = fullFit$opt$par
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=TRUE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=TRUE, reparameterizeVar=TRUE, varLogLambda=log(.000001))
# try removing the variance smoothness penalty: results look the same
initialParameters = fullFit$opt$par
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.1), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared = TRUE, estimateGpsShared = TRUE, includeGammaSpline = TRUE,
                      doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=FALSE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsVarGPS=1, nKnotsGamma=7, doSmoothnessPenalty=TRUE, conditionalGamma=TRUE, 
                      constrainMean=FALSE, maxCount=1, varSmoothnessPenalty=FALSE, nKnotsVar=5, paramInit=initialParameters, 
                      recompile=TRUE, reparameterizeVar=TRUE, varLogLambda=log(.000001))
# try increasing dStarGPS: x residuals look much better, but results turn non physical (mean turns skyhigh)
dStarGPS = 40000
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=11, nKnotsVarGPS=7, nKnotsGamma=7)
# increase penalty parameters, GPS knots: GPS mean doesn't need that many knots, x residuals biased shrunk to 0
dStarGPS = 40000
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(10),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=11, nKnotsVarGPS=11, nKnotsGamma=11)
# slightly decrease penalty parameters, reduce number of mean GPS knots: the mean becomes nonphysical again, y residuals get strange in the middle
dStarGPS = 40000
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(1),
                      sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, diffVar=TRUE, inflateVarLocking=TRUE, 
                      nKnotsMeanGPS=7, nKnotsVarGPS=11, nKnotsGamma=11)
testFit = fullFit
testFit$obj$par[length(testFit$obj$par)] = -5 # based on plotting the residuals
testFit$obj$theta = testFit$obj$par
testFit2 = fitModelTMB(testFit, fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                       G=G, debugPlotting=FALSE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                       sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE,
                       doMeanSpline=FALSE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                       dStarGPS=dStarGPS, diffMean=TRUE, inflateVarLocking=TRUE)
# when theshold is 40km and dStarGPS is 70000, gamma goes berserk, mean near 20, sd very small near 0.7:
# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE, 
#                       G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
#                       sharedSpatialProcess=FALSE, jointShared = FALSE, includeGammaSpline = FALSE, 
#                       doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar, 
#                       dStarGPS=dStarGPS, diffMean=TRUE)
# move threshold back to 25km, dStarGPS back to 25000, and make shared spatial process:
# x residuals seem unbiased and shrunk around 0, y seems underpredicted also with too low sds. mean and sd way too large
# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
#                       G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
#                       sharedSpatialProcess=TRUE, jointShared = FALSE, includeGammaSpline = FALSE,
#                       doMeanSpline=FALSE, doVarSpline=FALSE, diffGPSTaper=TRUE, dStar=dStar,
#                       dStarGPS=dStarGPS, diffMean=TRUE)
# most complex model without penalization: worked very well, mostly unbiased residuals, sd <= 2.6ish, mean shrinking to 0 moving southward, 
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=FALSE, includeGammaSpline=TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, constrainMean=TRUE)
testFit = fullFit
testFit$obj$par = testFit$minPar
testFit$obj$par[length(testFit$obj$par)] = -10
names(testFit$obj$par) = names(fullFit$obj$par)
testFit$obj$theta = testFit$obj$par
testFit2 = fitModelTMB(testFit, fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                       G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                       sharedSpatialProcess=TRUE, jointShared=FALSE, includeGammaSpline=TRUE,
                       doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                       dStarGPS=dStarGPS, diffMean=TRUE, constrainMean=TRUE)
# try the same thing but without constraining the mean: clearly not at the optimum based on gradient. residuals look awful and results are nonphysical
# [1] 474.4053
fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = FALSE,
                      G=G, debugPlotting=TRUE, logPenaltyPar=log(.0001), logDiffPenaltyPar=log(.0001),
                      sharedSpatialProcess=TRUE, jointShared=FALSE, includeGammaSpline=TRUE,
                      doMeanSpline=TRUE, doVarSpline=TRUE, diffGPSTaper=TRUE, dStar=dStar,
                      dStarGPS=dStarGPS, diffMean=TRUE, constrainMean=FALSE)