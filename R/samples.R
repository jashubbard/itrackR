epoch_samples <- function(obj,timevar,field='pa',epoch=c(-100,100))
{

  obj <- check_for_samples(obj) #if we haven't loaded sample data yet, this will do it

  # allepochs <- list()

  for(i in 1:length(obj$edfs))
  {

    #     justfile <- sub("^([^.]*).*", "\\1", basename(obj$edfs[[i]]))
    #     id <- as.numeric(gsub("([0-9]*).*","\\1",justfile))

    id <- edf2id(obj$edfs[[i]])

    header <- subset(obj$header,ID==id)
    events <- header[[timevar]][!is.na(header[[timevar]])]
    trials <- header$eyetrial[!is.na(header[[timevar]])]
    samples <- readRDS(obj$samples[[i]])

    thisepoch <- edfR::epoch.samples(events,as.data.frame(samples),sample.field=field,epoch=epoch,eyetrial=T)

    thisepoch$ID <- rep(id,length(trials))
    thisepoch$event <- timevar
    thisepoch$eyetrial <- trials
    obj$epochs$samples[[timevar]][[i]] <- thisepoch

    rm(samples)
    gc()

  }

  # obj$epochs[[timevar]] <- allepochs


  return(obj)

}

load_samples <- function(obj,outdir=NULL, force=F,parallel=TRUE,ncores = 2){


  allsamps <- list()
  samps <- data.table::data.table()
  alldata <- list()

  #check for directory to save .rds files into
  if(!is.null(outdir))
    obj$sample.dir <- outdir #if given as argument

  #otherwise create temporary directory
  else if(is.null(outdir) && (is.null(obj$sample.dir) || !dir.exists(obj$sample.dir))){
    tmpdir <- tempdir()
    dir.create(tmpdir,showWarnings = F)
    obj$sample.dir <- tmpdir

  }

  #default, run in parallel using foreach
  if(parallel){

    print('Loading samples (in parallel)...')


    #set up cluster with maximum number of cores
    ncores <- parallel::detectCores()
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    #run with dopar
    allsamps <- foreach::foreach(edf=iterators::iter(obj$edfs),.packages=c('edfR','data.table','itrackR')) %dopar%
    {
      fname <- load_sample_file(obj,edf)
      fname
    }

    parallel::stopCluster(cl)
  }
  else{

    #serial method - use lapply
    allsamps <- lapply(obj$edfs, function(x) load_sample_file(obj,x))

  }

  rm(alldata)
  rm(samps)



  obj$samples <- allsamps
  return(obj)

}

load_sample_file <- function(obj,edf){

   id <- edf2id(edf)
   fname <- file.path(obj$sample.dir,paste0(id,'_samp.rds'))

  if(file.exists(fname) && !force)
    samps <- readRDS(fname)
  else{
    samps <- edfR::edf.samples(edf,trials=T,eventmask=T)
    baseline <- obj$header$starttime[1]
    # samps <- timeshift(samps,baseline,'time')

    if(all(is.na(samps$paL))){
      samps[,(c("paL","gxL","gyL")) := NULL] #in-place delete of column
      data.table::setnames(samps,c("paR","gxR","gyR"),c("pa","gx","gy"))

    } else{

      samps[,(c("paR","gxR","gyR")) :=NULL]
      data.table::setnames(samps,c("paL","gxL","gyL"),c("pa","gx","gy"))

    }

    saveRDS(samps,fname,compress = T)
  }

  return(fname)
}



interpolate.blinks <- function(y,blinks)
{

  ## Locate blinks, blink starts and blink ends
  # A blink starts either when the blink variable goes from 0 to 1
  # or if its first value is 1. Similarly, a blink ends when the blink
  # variable goes from 1 to 0 or when the last value is 1.

  y.na <- is.na(y)
  if (all(!y.na)) return(y) # if no missing data, return y
  blink.start <- c(which.max(blinks), which(diff(blinks)==1) + 1)
  blink.start <- unique(blink.start) # remove eventual duplicates
  blink.end <- c(which(diff(blinks)==-1), max(which(blinks==1)))
  blink.end <- unique(blink.end) # remove eventual duplicates

  ## Interpolation
  n <- length(y)
  x.start <- pmax(blink.start-1,1)
  x.end <- pmin(blink.end+1,n)
  for (i in seq_along(x.start)) {
    xa <- x.start[i]
    xb <- x.end[i]
    if(!all(is.na(c(y[xa],y[xb]))))
      y[xa:xb] <- spline(x=c(xa,xb), y=c(y[xa],y[xb]), xout=xa:xb)$y
    else
      y[xa:xb] <- y[xa-1]
  }

  return(y)
}






remove_blinks <- function(obj, interpolate=FALSE)
{

  obj <- check_for_samples(obj)

  for(i in 1:length(obj$edfs))
  {

    samps<- readRDS(obj$samples[[i]])

    #code bad samples as "blinks"
    samps$blink[samps$gx==1e08 | samps$gy==1e08] <- 1

    blinks <- as.logical(samps$blink)
    samps$pa[blinks] <-NA
    samps$gx[blinks] <-NA
    samps$gy[blinks] <-NA

    if(interpolate){

      samps[, pa := interpolate.blinks(samps$pa,samps$blink)]
      samps[, gx := interpolate.blinks(samps$gx,samps$blink)]
      samps[, gy := interpolate.blinks(samps$gy,samps$blink)]
    }


    saveRDS(samps,obj$samples[[i]],compress = T)
    rm(samps)

  }

  return(obj)
}

get_all_epochs <- function(obj,epochname,baseline=NULL,baseline.method='percent',shape='wide',beh='all')
{

  window_start <- obj$epochs$samples[[epochname]][[1]]$epoch.window[1]
  window_end <- obj$epochs$samples[[epochname]][[1]]$epoch.window[2]

  epochs <- obj$epochs$samples[[epochname]]
  epochs <- lapply(epochs,function(x) cbind(x$ID,x$eyetrial,x$epochs))
  epochs <- do.call(rbind,epochs)

  if(!is.null(baseline))
  {
    epochs_b <- baseline_epochs(epochs[,-c(1,2)],baseline=baseline,method=baseline.method)
    epochs <- cbind(epochs[,1:2],epochs_b)
  }

  timenames <- paste0('t',seq(window_start,window_end-1))
  timenames <- gsub('-','_',timenames)
  colnames(epochs) <- c('ID','eyetrial',timenames)

  epochs <- as.data.frame(epochs)

  if(shape=='long'){

    epochs <- tidyr::gather_(epochs,'timepoint','value',timenames)
    epochs$timepoint <- gsub('t','',epochs$timepoint)
    epochs$timepoint <- as.numeric(gsub('_','-',epochs$timepoint))
  }

  if(!is.null(beh))
  {

    if(beh[1]=='all' || (is.logical(beh) && beh==TRUE))
      behnames <- names(obj$beh)
    else
      behnames <- unique(c('ID','eyetrial',beh))

    epochs <- dplyr::right_join(dplyr::select_(obj$beh,.dots=behnames),epochs,by=c('ID','eyetrial'))

    if(shape=='long')
      epochs <- dplyr::arrange(epochs,ID,eyetrial,timepoint)
    else
      epochs <- dplyr::arrange(epochs,ID,eyetrial)

  }

  return(epochs)
}


baseline_epochs <- function(epochs,baseline=c(1,100),method='percent'){


  mn <- rowMeans(epochs[,baseline[1]:baseline[2]],na.rm=T)

  epoch_b <- sweep(epochs,1,mn,FUN='-')
  epoch_b <- sweep(epoch_b,1,mn,FUN='/')

  return(epoch_b)
}




