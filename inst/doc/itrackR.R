## ----results="hide",message=FALSE----------------------------------------
library(itrackR)
datapath <- itrackR.data('path')
edf_files <- itrackR.data('edfs')
beh <- itrackR.data('beh')


## ------------------------------------------------------------------------
datapath


## ------------------------------------------------------------------------
edf_files

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(beh, 10))

## ------------------------------------------------------------------------
z <- itrackr(edfs=edf_files)

#Alternatively, we can provide the path and a search pattern to find all edfs in a certain folder:
#z <- itrackr(path=datapath, pattern='*.edf')

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(z$fixations, 10))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(z$saccades, 10))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(z$messages, 10))

## ------------------------------------------------------------------------
#generate coordinates for our ROIs
innercoords <- radialCoords(x=512,y=384,numpoints=6,radius=240);
outercoords <- radialCoords(512, 384,6, 280,starting_angle=30); #larger radius, starting w/ 30 degree offset

#specify rotations of ellipses
angles <- roiFlower(12)


## ---- fig.width=5, fig.height=4------------------------------------------

#make elliptical ROIs and plot them
z <- makeROIs(z,innercoords,shape='ellipse',xradius=60, yradius=120, angles=angles[c(1,3,5,7,9,11)])
plot.rois(z)


## ---- fig.width=5, fig.height=4------------------------------------------

#make elliptical ROIs and plot them
z <- makeROIs(z,outercoords,shape='ellipse',xradius=60, yradius=120, angles=angles[c(2,4,6,8,10,12)], names=7:12, append=T)
plot.rois(z)


## ----fig.width=5, fig.height=4-------------------------------------------

#coordinates have to be a matrix:
centercoords <- matrix(c(512,384),nrow=1)

z <- makeROIs(z,centercoords,shapes='circle',radius=65, names=13, append=T)
plot.rois(z)



## ---- fig.show='hold', fig.width=6, fig.height=4-------------------------
plot(z, zoom=TRUE)

## ------------------------------------------------------------------------
knitr::kable(head(z$header, 10))

## ------------------------------------------------------------------------
#find messages to use as our index variables (to merge with our behavioral data)
z<- set_index(z,varnames=c('Block','Trial'), patterns=c('^BLOCK [0-9]*','^TRIAL [0-9]*'), numeric.only=T)

#find messages that specify the onset of events, extract the timestamps
z <- find_messages(z,varnames=c('STIMONSET','RESPONSE'), patterns=c('STIMONSET','RESPONSE'), timestamp=T)
#merge with behavioral data
z <- add_behdata(z,beh,append=F)


## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(z$beh, 10))

## ------------------------------------------------------------------------
knitr::kable(head(z$header, 10))

## ----fig.width=6, fig.height=4-------------------------------------------
z <- drift_correct(z,vars='Block',threshold=15)
plot(z,zoom=T)

## ------------------------------------------------------------------------
z <- calcHits(z)

## ------------------------------------------------------------------------
knitr::kable(head(z$fixations, 10))

## ------------------------------------------------------------------------
z <- mapROIs(z,names=c('target','distractor'),indicators=c('Targetpos','Distractorpos'))

## ------------------------------------------------------------------------
knitr::kable(head(z$fixations, 10))

## ------------------------------------------------------------------------
fixes <- eyemerge(z,'fixations')

#including only some behavioral variables. ID and indexvars are always included
saccs <- eyemerge(z,'saccades',behdata=c('Task'))

#by default only mapped ROIs are included. Here we can include all 13 rois, plus the mapped ones
fixes_all <- eyemerge(z,'fixations',all.rois = T)


## ------------------------------------------------------------------------
#start at stimulus onset, going 700ms after that point. Bin the data into 25ms time bins.
z <- epoch_fixations(z,'target',event='STIMONSET',start = 0, end = 700, binwidth = 25)
z <- epoch_fixations(z,'distractor',event='STIMONSET',start = 0, end = 700, binwidth = 25)

## ----fig.width=7, fig.height=4-------------------------------------------
#plot the timeseries data for fixations to target, separately for the Conflict and Task conditions
plot.timeseries(z,event='STIMONSET',rois=c('target'),lines=c('Conflict'),cols='Task')

#plot fixations to target and disctractor for the same conditions
#you must specify 'roi' as one of the plotting variables for it to work.
plot.timeseries(z,event='STIMONSET',rois=c('target','distractor'),lines=c('roi','Conflict'),cols='Task')

#plot difference waves (target - distractor fixations). Plot on separate rows insetead of columns.
#This example doesn't make much sense because distractors aren't present on no-conflict trials
plot.timeseries(z,event='STIMONSET',rois=c('target','distractor'),lines=c('Conflict'),rows='Task',type='difference')


