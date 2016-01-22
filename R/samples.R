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

load_samples <- function(obj,outdir=NULL, force=F){



  allsamps <- list()
  samps <- data.table::data.table()
  alldata <- list()

  ncores <- parallel::detectCores()
  cl <- parallel::makeCluster((ncores/2)-1)
  doParallel::registerDoParallel(cl)


  if(!is.null(outdir))
    obj$sample.dir <- outdir
  else if(is.null(outdir) && is.null(obj$sample.dir))
    obj$sample.dir <- tempdir()



  # i <- 1
  allsamps <- foreach::foreach(edf=iterators::iter(obj$edfs),.packages=c('edfR','data.table','itrackR')) %dopar%
    # for(edf in obj$edfs)
  {
    # alldata <- edfR::edf.trials(edf,samples = T,eventmask = T)
    # print('Loading samples...')
    # samps <- data.table::data.table(alldata$samples)

    id <- edf2id(edf)
    fname <- file.path(obj$sample.dir,paste0(id,'_samp.rds'))

    if(file.exists(fname) && !force)
      samps <- readRDS(fname)
    else{
      samps <- data.table::as.data.table(edfR::edf.samples(edf,trials=T,eventmask=T))



    if(all(is.na(samps$paL))){
      samps[,(c("paL","gxL","gyL")) := NULL] #in-place delete of column
      data.table::setnames(samps,c("paR","gxR","gyR"),c("pa","gx","gy"))

      #       samps <- dplyr::select(samps,-paL,-gxL,gyL)
      #       samps <- dplyr::rename(samps,
      #                              pa = paR,
      #                              gx = gxR,
      #                              gy = gyR)
    } else{

      samps[,(c("paR","gxR","gyR")) :=NULL]
      data.table::setnames(samps,c("paL","gxL","gyL"),c("pa","gx","gy"))

      #       samps <- dplyr::select(samps,-paR,-gxR,gyR)
      #       samps <- dplyr::rename(samps,
      #                            pa = paL,
      #                            gx = gxL,
      #                            gy = gyL)
    }




    saveRDS(samps,fname,compress = T)
    }
    # allsamps[[i]] <- fname
    # rm(samps)
    # i <- i+1
    fname
  }

  rm(alldata)
  rm(samps)

  parallel::stopCluster(cl)

  obj$samples <- allsamps
  return(obj)

}

remove_blinks <- function(obj)
{

  obj <- check_for_samples(obj)

  for(i in 1:length(obj$edfs))
  {

    samps<- readRDS(obj$samples[[i]])

    blinks <- as.logical(samps$blink)
    samps$pa[blinks] <-NA
    samps$gx[blinks] <-NA
    samps$gy[blinks] <-NA

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




