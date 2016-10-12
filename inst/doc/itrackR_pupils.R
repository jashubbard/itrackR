## ----results="hide",message=FALSE----------------------------------------
library(itrackR)
datapath <- itrackr.data('path')
edf_files <- itrackr.data('edfs')
beh <- itrackr.data('beh')

z <- itrackr(edfs=edf_files)

## ------------------------------------------------------------------------
z <- find_messages(z,'STIMONSET','STIMONSET',timestamp = T)
z <- set_index(z,c('Block','Trial'),c('^BLOCK [0-9]*','TRIAL [0-9]*'),numeric.only = T)
z <- add_behdata(z,itrackr.data('beh'))

## ------------------------------------------------------------------------
z <- load_samples(z,outdir=datapath,parallel=T,ncores=2,force=T)

#see what your directory is..
z$sample.dir

## ------------------------------------------------------------------------

plot_samples(z,ID=104, nbins = 35, pages=2:3, events = T)

## ------------------------------------------------------------------------

z <- remove_blinks(z, interpolate=F)

## ------------------------------------------------------------------------

plot_samples(z,ID=104, nbins = 35, pages=2:3, events = T)

## ------------------------------------------------------------------------

#force it to reload
z <- load_samples(z,outdir=z$sample.dir,force=T)


## ------------------------------------------------------------------------
#the blinks are back!
plot_samples(z,ID=104, nbins = 35, pages=2:3, events = T)

## ------------------------------------------------------------------------
#the blinks are back!
z <- remove_blinks(z,interpolate=T)
plot_samples(z,ID=104, nbins = 35, pages=2:3, events = F, timestamp = 'STIMONSET')

## ------------------------------------------------------------------------

epochs <- get_sample_epochs(z,factors=c('Task','Conflict'), event='STIMONSET', epoch = c(-500,300), aggregate=T)

knitr::kable(head(epochs, 10))

## ------------------------------------------------------------------------
plot_sample_epochs(epochs,groups=c('Task'),colors = 'Conflict',rows='Task',aggregate=T)

## ------------------------------------------------------------------------

#get epochs again, this time baselining each trial 
epochs <- get_sample_epochs(z,factors=c('Task','Conflict'), event='STIMONSET', epoch = c(-500,300), aggregate=T, baseline = c(-500,-400))

#notice the plots are different
plot_sample_epochs(epochs,groups=c('Task'),colors = 'Conflict',rows='Task',aggregate=T)


