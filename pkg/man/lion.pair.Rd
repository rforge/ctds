\name{lion.pair}
\alias{lion.pair}
\docType{data}
\title{
Telemetry data for mountain lions.
}
\description{
Telemetry data for two mountain lions in the Front range cougar study.
An adult female, AF79, and her subadult cub, AM80.
}
\usage{data(lion.pair)}
\format{
  The format is:
List of 6
 $ lion1        :List of 7
  ..$ ID        : Factor w/ 57 levels "AF01","AF03",..: 31
  ..$ Sex       : Factor w/ 2 levels "Female","Male": 1
  ..$ locs.orig :'data.frame':	540 obs. of  3 variables:
  .. ..$ UTM_X   : num [1:540] 458189 458212 457834 457840 458567 ...
  .. ..$ UTM_Y   : num [1:540] 4384107 4384188 4384561 4384554 4384293 ...
  .. ..$ datetime: num [1:540] 4095 4095 4095 4096 4096 ...
  ..$ locs      :'data.frame':	98 obs. of  3 variables:
  .. ..$ UTM_X   : num [1:98] 458189 458212 457834 457840 458567 ...
  .. ..$ UTM_Y   : num [1:98] 4384107 4384188 4384561 4384554 4384293 ...
  .. ..$ datetime: num [1:98] 4095 4095 4095 4096 4096 ...
  ..$ crawl.fit :List of 34
  .. ..$ par          : num [1:4] 0.916 0.916 7.482 15.126
  .. ..$ estPar       : num [1:2] 7.48 15.13
  .. ..$ se           : num [1:4] NA NA 0.0508 297.3302
  .. ..$ ci           : num [1:4, 1:2] NA NA 7.38 -567.64 NA ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr [1:2] "L" "U"
  .. ..$ Cmat         : num [1:4, 1:4] NA NA NA NA NA NA NA NA NA NA ...
  .. ..$ loglik       : num -1386
  .. ..$ aic          : num 2777
  .. ..$ initial.state:List of 4
  .. .. ..$ a1.x: num [1:2] 458189 0
  .. .. ..$ a1.y: num [1:2] 4384107 0
  .. .. ..$ P1.x: num [1:2, 1:2] 1 0 0 1
  .. .. ..$ P1.y: num [1:2, 1:2] 1 0 0 1
  .. ..$ coord        : chr [1:2] "UTM_X" "UTM_Y"
  .. ..$ fixPar       : num [1:4] 0.916 0.916 NA NA
  .. ..$ convergence  : int 0
  .. ..$ message      : chr "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  .. ..$ stop.model   : NULL
  .. ..$ random.drift : logi FALSE
  .. ..$ mov.model    :Class 'formula' length 2 ~1
  .. .. .. ..- attr(*, ".Environment")=<environment: 0xaef3960> 
  .. ..$ err.model    :List of 2
  .. .. ..$ x:Class 'formula' length 2 ~1
  .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. .. ..$ y:Class 'formula' length 2 ~1
  .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. ..$ n.par        : int 4
  .. ..$ nms          : chr [1:4] "ln tau.x (Intercept)" "ln tau.y (Intercept)" "ln sigma (Intercept)" "ln beta (Intercept)"
  .. ..$ n.mov        : int 1
  .. ..$ n.errX       : int 1
  .. ..$ n.errY       : int 1
  .. ..$ mov.mf       : num [1:98, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. .. ..- attr(*, "assign")= int 0
  .. ..$ err.mfX      : num [1:98, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. ..$ err.mfY      : num [1:98, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. ..$ stop.mf      : NULL
  .. ..$ polar.coord  : logi FALSE
  .. ..$ Time.name    : chr "datetime"
  .. ..$ init         :List of 5
  .. .. ..$ par        : num [1:2] 7.64 3.44
  .. .. ..$ value      : num 2778
  .. .. ..$ counts     : Named int [1:2] 50 NA
  .. .. .. ..- attr(*, "names")= chr [1:2] "function" "gradient"
  .. .. ..$ convergence: int 0
  .. .. ..$ message    : NULL
  .. ..$ data         :'data.frame':	98 obs. of  3 variables:
  .. .. ..$ UTM_X   : num [1:98] 458189 458212 457834 457840 458567 ...
  .. .. ..$ UTM_Y   : num [1:98] 4384107 4384188 4384561 4384554 4384293 ...
  .. .. ..$ datetime: num [1:98] 4095 4095 4095 4096 4096 ...
  .. ..$ lower        : num [1:2] -Inf -Inf
  .. ..$ upper        : num [1:2] Inf Inf
  .. ..$ prior        : NULL
  .. ..$ need.hess    : num 1
  .. ..$ runTime      :Class 'difftime'  atomic [1:1] 0.754
  .. .. .. ..- attr(*, "units")= chr "secs"
  .. ..- attr(*, "class")= chr "crwFit"
  ..$ pred.times: num [1:1015] 4095 4095 4095 4095 4095 ...
  ..$ pred.obj  :List of 9
  .. ..$ originalData:'data.frame':	1112 obs. of  4 variables:
  .. .. ..$ datetime: num [1:1112] 4095 4095 4095 4095 4095 ...
  .. .. ..$ locType : chr [1:1112] "p" "p" "p" "p" ...
  .. .. ..$ UTM_X   : num [1:1112] 458189 NA NA NA NA ...
  .. .. ..$ UTM_Y   : num [1:1112] 4384107 NA NA NA NA ...
  .. ..$ alpha.hat.y :'data.frame':	1112 obs. of  2 variables:
  .. .. ..$ mu.y: num [1:1112] 4384107 4384120 4384134 4384147 4384161 ...
  .. .. ..$ nu.y: num [1:1112] 5.58e-11 6.53e+02 6.53e+02 6.53e+02 6.53e+02 ...
  .. ..$ alpha.hat.x :'data.frame':	1112 obs. of  2 variables:
  .. .. ..$ mu.x: num [1:1112] 458189 458193 458197 458201 458204 ...
  .. .. ..$ nu.x: num [1:1112] 1.54e-11 1.80e+02 1.80e+02 1.80e+02 1.80e+02 ...
  .. ..$ V.hat.y     : num [1:2, 1:2, 1:1112] 1 0 0 1 54734 ...
  .. ..$ V.hat.x     : num [1:2, 1:2, 1:1112] 1 0 0 1 54734 ...
  .. ..$ speed       :'data.frame':	1112 obs. of  2 variables:
  .. .. ..$ ln.speed    : num [1:1112] -23.57 6.52 6.52 6.52 6.52 ...
  .. .. ..$ var.ln.speed: num [1:1112] 2.98e+20 1.28e+07 1.28e+07 1.28e+07 1.28e+07 ...
  .. ..$ loglik      : num -1386
  .. ..$ useAvail    :List of 2
  .. .. ..$ use  : num [1:1014, 1:5] 4095 4095 4095 4095 4095 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:5] "datetime" "meanUse.x" "meanUse.y" "varUse.x" ...
  .. .. ..$ avail: num [1:1014, 1:5] 4095 4095 4095 4095 4095 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:5] "datetime" "meanAvail.x" "meanAvail.y" "varAvail.x" ...
  .. ..$ fit.test    :'data.frame':	1112 obs. of  8 variables:
  .. .. ..$ shock.x      : num [1:1112] -5.7e-05 NA NA NA NA ...
  .. .. ..$ predObs.x    : num [1:1112] 458189 NA NA NA NA ...
  .. .. ..$ Vshock.x     : num [1:1112] 7.25 NA NA NA NA ...
  .. .. ..$ shock.y      : num [1:1112] -0.000207 NA NA NA NA ...
  .. .. ..$ predObs.y    : num [1:1112] 4384107 NA NA NA NA ...
  .. .. ..$ Vshock.y     : num [1:1112] 7.25 NA NA NA NA ...
  .. .. ..$ outlier.chisq: num [1:1112] 6.36e-09 NaN NaN NaN NaN ...
  .. .. ..$ naive.p.val  : num [1:1112] 1 NaN NaN NaN NaN ...
  .. ..- attr(*, "flat")= logi FALSE
  .. ..- attr(*, "coord")= Named chr [1:2] "UTM_X" "UTM_Y"
  .. .. ..- attr(*, "names")= chr [1:2] "x" "y"
  .. ..- attr(*, "random.drift")= logi FALSE
  .. ..- attr(*, "polar.coord")= logi FALSE
  .. ..- attr(*, "Time.name")= chr "datetime"
  .. ..- attr(*, "class")= chr [1:2] "list" "crwPredict"
 $ lion2        :List of 8
  ..$ ID        : Factor w/ 57 levels "AF01","AF03",..: 54
  ..$ Sex       : Factor w/ 2 levels "Female","Male": 2
  ..$ locs.orig :'data.frame':	885 obs. of  3 variables:
  .. ..$ UTM_X   : num [1:885] 458318 458613 458024 458024 458023 ...
  .. ..$ UTM_Y   : num [1:885] 4384122 4384338 4384590 4384604 4384593 ...
  .. ..$ datetime: num [1:885] 4095 4095 4095 4096 4096 ...
  ..$ locs      :'data.frame':	140 obs. of  3 variables:
  .. ..$ UTM_X   : num [1:140] 458318 458613 458024 458024 458023 ...
  .. ..$ UTM_Y   : num [1:140] 4384122 4384338 4384590 4384604 4384593 ...
  .. ..$ datetime: num [1:140] 4095 4095 4095 4096 4096 ...
  ..$ crawl.fit :List of 34
  .. ..$ par          : num [1:4] 0.916 0.916 7.026 3.81
  .. ..$ estPar       : num [1:2] 7.03 3.81
  .. ..$ se           : num [1:4] NA NA 0.0659 0.499
  .. ..$ ci           : num [1:4, 1:2] NA NA 6.9 2.83 NA ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr [1:2] "L" "U"
  .. ..$ Cmat         : num [1:4, 1:4] NA NA NA NA NA NA NA NA NA NA ...
  .. ..$ loglik       : num -1803
  .. ..$ aic          : num 3610
  .. ..$ initial.state:List of 4
  .. .. ..$ a1.x: num [1:2] 458318 0
  .. .. ..$ a1.y: num [1:2] 4384122 0
  .. .. ..$ P1.x: num [1:2, 1:2] 1 0 0 1
  .. .. ..$ P1.y: num [1:2, 1:2] 1 0 0 1
  .. ..$ coord        : chr [1:2] "UTM_X" "UTM_Y"
  .. ..$ fixPar       : num [1:4] 0.916 0.916 NA NA
  .. ..$ convergence  : int 0
  .. ..$ message      : chr "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  .. ..$ stop.model   : NULL
  .. ..$ random.drift : logi FALSE
  .. ..$ mov.model    :Class 'formula' length 2 ~1
  .. .. .. ..- attr(*, ".Environment")=<environment: 0xbf57940> 
  .. ..$ err.model    :List of 2
  .. .. ..$ x:Class 'formula' length 2 ~1
  .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. .. ..$ y:Class 'formula' length 2 ~1
  .. .. .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
  .. ..$ n.par        : int 4
  .. ..$ nms          : chr [1:4] "ln tau.x (Intercept)" "ln tau.y (Intercept)" "ln sigma (Intercept)" "ln beta (Intercept)"
  .. ..$ n.mov        : int 1
  .. ..$ n.errX       : int 1
  .. ..$ n.errY       : int 1
  .. ..$ mov.mf       : num [1:140, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. .. ..- attr(*, "assign")= int 0
  .. ..$ err.mfX      : num [1:140, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. ..$ err.mfY      : num [1:140, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr "(Intercept)"
  .. ..$ stop.mf      : NULL
  .. ..$ polar.coord  : logi FALSE
  .. ..$ Time.name    : chr "datetime"
  .. ..$ init         :List of 5
  .. .. ..$ par        : num [1:2] 7 3.94
  .. .. ..$ value      : num 3606
  .. .. ..$ counts     : Named int [1:2] 200 NA
  .. .. .. ..- attr(*, "names")= chr [1:2] "function" "gradient"
  .. .. ..$ convergence: int 0
  .. .. ..$ message    : NULL
  .. ..$ data         :'data.frame':	140 obs. of  3 variables:
  .. .. ..$ UTM_X   : num [1:140] 458318 458613 458024 458024 458023 ...
  .. .. ..$ UTM_Y   : num [1:140] 4384122 4384338 4384590 4384604 4384593 ...
  .. .. ..$ datetime: num [1:140] 4095 4095 4095 4096 4096 ...
  .. ..$ lower        : num [1:2] -Inf -Inf
  .. ..$ upper        : num [1:2] Inf Inf
  .. ..$ prior        : NULL
  .. ..$ need.hess    : num 1
  .. ..$ runTime      :Class 'difftime'  atomic [1:1] 0.855
  .. .. .. ..- attr(*, "units")= chr "secs"
  .. ..- attr(*, "class")= chr "crwFit"
  ..$ pred.times: num [1:1015] 4095 4095 4095 4095 4095 ...
  ..$ pred.obj  :List of 9
  .. ..$ originalData:'data.frame':	1154 obs. of  4 variables:
  .. .. ..$ datetime: num [1:1154] 4095 4095 4095 4095 4095 ...
  .. .. ..$ locType : chr [1:1154] "p" "p" "p" "p" ...
  .. .. ..$ UTM_X   : num [1:1154] 458318 NA NA NA NA ...
  .. .. ..$ UTM_Y   : num [1:1154] 4384122 NA NA NA NA ...
  .. ..$ alpha.hat.y :'data.frame':	1154 obs. of  2 variables:
  .. .. ..$ mu.y: num [1:1154] 4384122 4384137 4384170 4384208 4384250 ...
  .. .. ..$ nu.y: num [1:1154] 3.61e-05 1.26e+03 1.75e+03 1.95e+03 2.03e+03 ...
  .. ..$ alpha.hat.x :'data.frame':	1154 obs. of  2 variables:
  .. .. ..$ mu.x: num [1:1154] 458318 458346 458406 458476 458544 ...
  .. .. ..$ nu.x: num [1:1154] 6.79e-05 2.35e+03 3.21e+03 3.40e+03 3.09e+03 ...
  .. ..$ V.hat.y     : num [1:2, 1:2, 1:1154] 1 0 0 1 3341 ...
  .. ..$ V.hat.x     : num [1:2, 1:2, 1:1154] 1 0 0 1 3341 ...
  .. ..$ speed       :'data.frame':	1154 obs. of  2 variables:
  .. .. ..$ ln.speed    : num [1:1154] -9.47 7.89 8.21 8.27 8.21 ...
  .. .. ..$ var.ln.speed: num [1:1154] 1.69e+08 2.69 1.36 1.11 1.28 ...
  .. ..$ loglik      : num -1803
  .. ..$ useAvail    :List of 2
  .. .. ..$ use  : num [1:1014, 1:5] 4095 4095 4095 4095 4095 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:5] "datetime" "meanUse.x" "meanUse.y" "varUse.x" ...
  .. .. ..$ avail: num [1:1014, 1:5] 4095 4095 4095 4095 4095 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:5] "datetime" "meanAvail.x" "meanAvail.y" "varAvail.x" ...
  .. ..$ fit.test    :'data.frame':	1154 obs. of  8 variables:
  .. .. ..$ shock.x      : num [1:1154] -0.00309 NA NA NA NA ...
  .. .. ..$ predObs.x    : num [1:1154] 458318 NA NA NA NA ...
  .. .. ..$ Vshock.x     : num [1:1154] 7.25 NA NA NA NA ...
  .. .. ..$ shock.y      : num [1:1154] -0.00163 NA NA NA NA ...
  .. .. ..$ predObs.y    : num [1:1154] 4384122 NA NA NA NA ...
  .. .. ..$ Vshock.y     : num [1:1154] 7.25 NA NA NA NA ...
  .. .. ..$ outlier.chisq: num [1:1154] 1.69e-06 NaN NaN NaN NaN ...
  .. .. ..$ naive.p.val  : num [1:1154] 1 NaN NaN NaN NaN ...
  .. ..- attr(*, "flat")= logi FALSE
  .. ..- attr(*, "coord")= Named chr [1:2] "UTM_X" "UTM_Y"
  .. .. ..- attr(*, "names")= chr [1:2] "x" "y"
  .. ..- attr(*, "random.drift")= logi FALSE
  .. ..- attr(*, "polar.coord")= logi FALSE
  .. ..- attr(*, "Time.name")= chr "datetime"
  .. ..- attr(*, "class")= chr [1:2] "list" "crwPredict"
  ..$ sim.obj   :List of 32
  .. ..$ x            : num [1:1154] 458318 9999 9999 9999 9999 ...
  .. ..$ y            : num [1:1154] 4384122 9999 9999 9999 9999 ...
  .. ..$ loctype      : num [1:1154] 0 1 1 1 1 1 1 0 1 1 ...
  .. ..$ P1.y         : num [1:2, 1:2] 1 0 0 1
  .. ..$ P1.x         : num [1:2, 1:2] 1 0 0 1
  .. ..$ a1.x         : num [1:2] 458318 0
  .. ..$ a1.y         : num [1:2] 4384122 0
  .. ..$ n.errX       : int 1
  .. ..$ n.errY       : int 1
  .. ..$ n.mov        : int 1
  .. ..$ delta        : num [1:1154] 0.0208 0.0208 0.0208 0.0208 0.0208 ...
  .. ..$ lonAdj       : num [1:1154] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ driftMod     : logi FALSE
  .. ..$ stopMod      : logi FALSE
  .. ..$ stop.mf      : NULL
  .. ..$ err.mfX      : num [1:1154, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ err.mfY      : num [1:1154, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ mov.mf       : num [1:1154, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ fixPar       : num [1:4] 0.916 0.916 NA NA
  .. ..$ Cmat         : num [1:4, 1:4] NA NA NA NA NA NA NA NA NA NA ...
  .. ..$ predType     : chr [1:1154] "p" "p" "p" "p" ...
  .. ..$ par          : num [1:4] 0.916 0.916 7.026 3.81
  .. ..$ nms          : chr [1:4] "ln tau.x (Intercept)" "ln tau.y (Intercept)" "ln sigma (Intercept)" "ln beta (Intercept)"
  .. ..$ N            : int 1154
  .. ..$ lower        : num [1:2] -Inf -Inf
  .. ..$ upper        : num [1:2] Inf Inf
  .. ..$ loglik       : num -1803
  .. ..$ datetime     : num [1:1154] 4095 4095 4095 4095 4095 ...
  .. ..$ coord        : chr [1:2] "UTM_X" "UTM_Y"
  .. ..$ Time.name    : chr "datetime"
  .. ..$ prior        : NULL
  .. ..$ thetaSampList:List of 1
  .. .. ..$ : num [1:1000, 1:7] 0.000944 0.000866 0.000776 0.000938 0.000912 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:7] "w" "lik" "prop.lik" "ln tau.x (Intercept)" ...
  .. .. .. ..- attr(*, "effSamp")= num 137
  .. .. .. ..- attr(*, "method")= chr "IS"
  .. .. .. ..- attr(*, "numLikEval")= num 1000
  .. ..- attr(*, "class")= chr "crwSimulator"
 $ mintime      : num 4096
 $ maxtime      : num 4114
 $ extent       :Formal class 'Extent' [package "raster"] with 4 slots
  .. ..@ xmin: num 398400
  .. ..@ xmax: num 469197
  .. ..@ ymin: num 4355586
  .. ..@ ymax: num 4438958
 $ stack.cropped:Formal class 'RasterStack' [package "raster"] with 11 slots
  .. ..@ filename  : chr ""
  .. ..@ layers    :List of 7
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi FALSE
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "city"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi NA
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "ag"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi NA
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "grass"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi NA
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "shrub"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi NA
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "forest"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi NA
  .. .. .. .. .. .. ..@ attributes:List of 1
  .. .. .. .. .. .. .. ..$ : NULL
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 0
  .. .. .. .. .. .. ..@ max       : num 1
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "bare"
  .. .. .. .. ..@ z         : list()
  .. .. ..$ :Formal class 'RasterLayer' [package "raster"] with 13 slots
  .. .. .. .. ..@ file      :Formal class '.RasterFile' [package "raster"] with 11 slots
  .. .. .. .. .. .. ..@ name        : chr ""
  .. .. .. .. .. .. ..@ datanotation: chr "FLT4S"
  .. .. .. .. .. .. ..@ byteorder   : chr "little"
  .. .. .. .. .. .. ..@ nodatavalue : num -Inf
  .. .. .. .. .. .. ..@ nbands      : int 1
  .. .. .. .. .. .. ..@ bandorder   : chr "BIL"
  .. .. .. .. .. .. ..@ offset      : int 0
  .. .. .. .. .. .. ..@ toptobottom : logi TRUE
  .. .. .. .. .. .. ..@ blockrows   : int 0
  .. .. .. .. .. .. ..@ blockcols   : int 0
  .. .. .. .. .. .. ..@ driver      : chr ""
  .. .. .. .. ..@ data      :Formal class '.SingleLayerData' [package "raster"] with 12 slots
  .. .. .. .. .. .. ..@ values    : num [1:15260] 2997 2976 2958 2929 2901 ...
  .. .. .. .. .. .. ..@ offset    : num 0
  .. .. .. .. .. .. ..@ gain      : num 1
  .. .. .. .. .. .. ..@ inmemory  : logi TRUE
  .. .. .. .. .. .. ..@ fromdisk  : logi FALSE
  .. .. .. .. .. .. ..@ isfactor  : logi FALSE
  .. .. .. .. .. .. ..@ attributes: list()
  .. .. .. .. .. .. ..@ haveminmax: logi TRUE
  .. .. .. .. .. .. ..@ min       : num 2263
  .. .. .. .. .. .. ..@ max       : num 3597
  .. .. .. .. .. .. ..@ band      : int 1
  .. .. .. .. .. .. ..@ unit      : chr ""
  .. .. .. .. ..@ legend    :Formal class '.RasterLegend' [package "raster"] with 5 slots
  .. .. .. .. .. .. ..@ type      : chr(0) 
  .. .. .. .. .. .. ..@ values    : logi(0) 
  .. .. .. .. .. .. ..@ color     : logi(0) 
  .. .. .. .. .. .. ..@ names     : logi(0) 
  .. .. .. .. .. .. ..@ colortable: logi(0) 
  .. .. .. .. ..@ history   : list()
  .. .. .. .. ..@ title     : chr(0) 
  .. .. .. .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. .. .. .. ..@ xmin: num 451524
  .. .. .. .. .. .. ..@ xmax: num 465524
  .. .. .. .. .. .. ..@ ymin: num 4378808
  .. .. .. .. .. .. ..@ ymax: num 4389708
  .. .. .. .. ..@ rotated   : logi FALSE
  .. .. .. .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. .. .. .. ..@ geotrans: num(0) 
  .. .. .. .. .. .. ..@ transfun:function ()  
  .. .. .. .. ..@ ncols     : int 140
  .. .. .. .. ..@ nrows     : int 109
  .. .. .. .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. .. .. ..@ projargs: chr " +proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. .. .. .. ..@ layernames: chr "r_1"
  .. .. .. .. ..@ z         : list()
  .. ..@ title     : chr(0) 
  .. ..@ extent    :Formal class 'Extent' [package "raster"] with 4 slots
  .. .. .. ..@ xmin: num 451524
  .. .. .. ..@ xmax: num 465524
  .. .. .. ..@ ymin: num 4378808
  .. .. .. ..@ ymax: num 4389708
  .. ..@ rotated   : logi FALSE
  .. ..@ rotation  :Formal class '.Rotation' [package "raster"] with 2 slots
  .. .. .. ..@ geotrans: num(0) 
  .. .. .. ..@ transfun:function ()  
  .. ..@ ncols     : int 140
  .. ..@ nrows     : int 109
  .. ..@ crs       :Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. ..@ projargs: chr "+proj=utm +zone=13 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"
  .. ..@ layernames: chr [1:7] "city" "ag" "grass" "shrub" ...
  .. ..@ z         : list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(lion.pair)
## maybe str(lion.pair) ; plot(lion.pair) ...
}
\keyword{datasets}
