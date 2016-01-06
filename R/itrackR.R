

#class constructor
itrackr <- function(edfs = NULL,path=NULL,pattern=NULL,resolution=c(1024,768),datadir=tempdir())
{

  me <- list(
              edfs = edfs,
              ids = NULL,
              idvar = 'ID',
              indexvars = NULL,
              timevars = NULL,
              resolution = resolution,
              header = data.frame,
              samples = list(),
              fixations = data.frame,
              saccades = data.frame,
              blinks = data.frame,
              messages = data.frame,
              epochs = list(),
              beh = data.frame
  )


  class(me) <- append('itrackR',class(me))

  if(!is.null(me$edfs))
    {
      path <- edfs
      pattern <- NULL
    }

  if(!is.null(path))
    me <- load_edfs(me,path=path,pattern=pattern)


  return(me)

}


load_edfs <- function(obj,path='.',pattern='*.edf',recursive = FALSE){

  fields <- c('fixations','saccades','blinks','messages','header')


  allbatch <- edfR::edf.batch(path,pattern=pattern,samples=F,do.plot=F)

  alldata <- edfR::combine.eyedata(allbatch,fields=fields)

  obj$header <- alldata$header
  obj$fixations <- alldata$fixations
  obj$fixations$fixation_key <- 1:nrow(obj$fixations)
  obj$saccades <- alldata$saccades
  obj$saccades$saccade_key <- 1:nrow(obj$saccades)
  obj$blinks <-alldata$blinks
  obj$blinks$blink_key <- 1:nrow(obj$blinks)
  obj$messages <- alldata$messages
  obj$messages$message_key <- 1:nrow(obj$messages)

#   if(samples)
#     {
#       allsamples <- load_samples(obj)
#       obj$samples <- allsamples$samples
#     }

  if(is.null(obj$edfs))
    obj$edfs <- unlist(lapply(tmp,function(x) x$filename))

  obj$edfs <- unlist(lapply(obj$edfs,make_fullpath))


  return(obj)

}

check_for_samples <- function(obj)
{

  if(length(obj$samples)<length(obj$edfs))
  {
    print('Sample data not loaded yet. Loading now (this make take a while)...')
    obj <- load_samples(obj)
  }



 return(obj)

}


make_fullpath <- function(fname){

  if(dirname(fname)=='.')
  {
    pathdir <- getwd()
    fullname <- file.path(pathdir,fname)
  }
  else
    fullname <- path.expand(fname)

  return(fullname)
}


plot <- function(obj,zoom=FALSE,crosshairs=TRUE,rois=TRUE,which='all',names=FALSE){
  UseMethod('plot')

}

plot.itrackR <- function(obj,zoom=FALSE,crosshairs=TRUE,rois=TRUE,which='all',names=FALSE){


  if(rois && !is.null(obj$rois)){
    df <- rois2df(obj)

    if(which !='all')
      df <- subset(df,name %in% which)


    p <- ggplot2::ggplot() + ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)+
      ggplot2::geom_point(data=obj$fixations,ggplot2::aes(x=gavx,y=gavy),color='yellow',size=0.7)

    if(names)
      p <- p + ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),aes(x=xcenter,y=ycenter,label=name),color='white')
  }
  else
    p <- ggplot2::ggplot(obj$fixations,ggplot2::aes(x=gavx,y=gavy)) + ggplot2::geom_point(color='yellow',size=0.7)


  if(zoom)
    p <- p+ggplot2::xlim(xlim=c(0,obj$resolution[1])) + ggplot2::ylim(ylim=c(obj$resolution[2],0))
  else
    p <- p + ggplot2::xlim(xlim=c(0,max(obj$fixations$gavx))) + ggplot2::ylim(ylim=c(max(obj$fixations$gavy),0))

    p <- p + ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank())


    if(crosshairs){
      p <- p +
        ggplot2::geom_hline(yintercept=round(z$resolution[2]/2),color='red',size=0.3,linetype='dashed') +
        ggplot2::geom_vline(xintercept=round(z$resolution[1]/2),color='red',size=0.3,linetype='dashed')
    }


  p <- p + ggplot2::facet_wrap(~ID)


  p
  return(p)
}


set_index <- function(obj,varnames,patterns,numeric.only=FALSE)
{
  obj <- find_messages(obj,varnames,patterns,numeric.only)

  obj$indexvars <- varnames

  return(obj)
}



find_messages <- function(obj,varnames,patterns,numeric.only = FALSE,timestamp=FALSE)
  {

  msg_index <- dplyr::distinct(dplyr::select(obj$header,ID,eyetrial))

  for(i in 1:length(varnames)){

    find_msg <- grepl(patterns[i],obj$messages$message)

    tmp <- obj$messages[find_msg,c('ID','eyetrial')]

    if(numeric.only)
      msgs <- as.numeric(gsub("[^0-9]", "",obj$messages$message[find_msg]))
    else if(!numeric.only && timestamp)
      msgs <- obj$messages$sttime[find_msg]
    else
      msgs <- obj$messages$message[find_msg]

    tmp[varnames[i]] <-msgs

    msg_index <- dplyr::left_join(msg_index,tmp,by=c('ID','eyetrial'))

  }


  matches <- names(obj$header) %in% varnames
  obj$header <- obj$header[!matches]


  obj$header <- dplyr::left_join(obj$header,msg_index,by=c('ID','eyetrial'))

  if(timestamp)
    obj$timevars <- varnames

  return(obj)

}


add_behdata <- function(obj,beh,append=FALSE){

  if(append){

    obj$beh <- merge(obj$beh,beh,by=c('ID',obj$indexvars))

    }
  else{

    eyedata <- obj$header[c('ID','eyetrial',obj$indexvars)]
    eyedata <- cbind(eyedata,obj$header[obj$timevars] - obj$header$starttime)

    behmerged <- merge(beh,eyedata,by=c('ID',obj$indexvars),all.x=T)
    obj$beh <- behmerged

    }

  return(obj)
}


epoch_samples <- function(obj,timevar,field='pa',epoch=c(-100,100),cleanup=F)
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
    obj$epochs[[timevar]][[i]] <- thisepoch

    rm(samples)
    gc()

  }

  # obj$epochs[[timevar]] <- allepochs


  return(obj)

}


get_subdata <- function(obj,id,fields = c('header','samples'))
{
  subobj <- itrackr()

  for(f in fields)
  {

    subobj[[f]] <- subset(obj[[f]],ID==id)

  }

  return(subobj)




}


load_samples <- function(obj,outdir=tempdir()){

  allsamps <- list()
  samps <- data.table::data.table()
  alldata <- list()

  ncores <- parallel::detectCores()
  cl <- parallel::makeCluster((ncores/2)-1)
  doParallel::registerDoParallel(cl)

  # i <- 1
  allsamps <- foreach::foreach(edf=iterators::iter(obj$edfs),.packages=c('edfR','data.table','itrackR')) %dopar%
  # for(edf in obj$edfs)
  {
    # alldata <- edfR::edf.trials(edf,samples = T,eventmask = T)

    # samps <- data.table::data.table(alldata$samples)
      samps <- data.table::as.data.table(edfR::edf.samples(edf,trials=T,eventmask=T))

    if(all(is.na(samps$paL))){
      samps[,c("paL","gxL","gyL") := NULL] #in-place delete of column
      data.table::setnames(samps,c("paR","gxR","gyR"),c("pa","gx","gy"))

#       samps <- dplyr::select(samps,-paL,-gxL,gyL)
#       samps <- dplyr::rename(samps,
#                              pa = paR,
#                              gx = gxR,
#                              gy = gyR)
      }
    else{

      samps[,c("paR","gxR","gyR") :=NULL]
      data.table::setnames(samps,c("paL","gxL","gyL"),c("pa","gx","gy"))

#       samps <- dplyr::select(samps,-paR,-gxR,gyR)
#       samps <- dplyr::rename(samps,
#                            pa = paL,
#                            gx = gxL,
#                            gy = gyL)
    }


    id <- edf2id(edf)
    fname <- file.path(outdir,paste0(id,'_samp.rds'))
    saveRDS(samps,fname,compress = T)
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

get_all_epochs <- function(obj,epochname,baseline=NULL,baseline.method='percent',shape='wide',beh=NULL)
{

  window_start <- obj$epochs[[epochname]][[1]]$epoch.window[1]
  window_end <- obj$epochs[[epochname]][[1]]$epoch.window[2]

  epochs <- obj$epochs[[epochname]]
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

    if(is.logical(beh) && beh==TRUE)
      behnames <- names(obj$beh)
    else
      behnames <- unique(c('ID','eyetrial',beh))

    epochs <- dplyr::left_join(epochs,dplyr::select_(obj$beh,.dots=behnames),by=c('ID','eyetrial'))

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

plot_random_epochs <- function(epochs,n=100)
{

  if('timepoint' %in% names(epochs))
  {
    epochs <- dplyr::select(epochs,ID,eyetrial,timepoint,value)
    epochs <- tidyr::spread(epochs,timepoint,value)
  }


  epochs <- dplyr::select(epochs,matches('^t[_0-9]'))

  tmp <- dplyr::sample_n(epochs,n)
  matplot(t(as.matrix(tmp)),type='l',lty=1)

}

remove_ext <- function(filename){
  justfile <- sub("^([^.]*).*", "\\1", basename(filename))
  return(justfile)
}

edf2id <- function(edf){

  justfile <- remove_ext(edf)
  id <- as.numeric(gsub("([0-9]*).*","\\1",justfile))
  return(id)
}



makeROIs <- function(obj,coords,shapes='circle',radius=0,xradius=0,yradius=0,angles=NULL,names=NULL,append=F){

  if(length(shapes)==1 && nrow(coords)>1)
    shapes <- rep(shapes[1],nrow(coords))

  if(length(radius)==1 && nrow(coords)>1)
    radius <- rep(radius[1],nrow(coords))

  if(length(xradius)==1 && nrow(coords)>1)
    xradius <- rep(xradius[1],nrow(coords))

  if(length(yradius)==1 && nrow(coords)>1)
    yradius <- rep(yradius[1],nrow(coords))



  if(length(obj$rois)==0 || append==F){
    start_roi=1
    obj$rois <- list()
  }
  else if(length(obj$rois)>0 && append==T){
   start_roi = length(obj$rois)+1
  }


  if(is.null(names))
    names <- start_roi:nrow(coords)

  roipos <- start_roi


 for(i in 1:nrow(coords)){

    tmpROI <- list()

    if(shapes[[i]]=='circle')
      tmpROI$roi <- disc(radius=radius[[i]],centre=coords[i,])
    else if(shapes[[i]]=='ellipse')
      tmpROI$roi <- ellipse(xradius[[i]],yradius[[i]],centre=coords[i,],phi=angles[[i]])

    tmpROI$name <- names[[i]]
    tmpROI$shape <- shapes[[i]]
    tmpROI$center <- coords[i,]
    tmpROI$radius <- radius[[i]]
    tmpROI$xradius <- xradius[[i]]
    tmpROI$yradius <- yradius[[i]]

    obj$rois[[roipos]] <- tmpROI
    roipos <- roipos + 1

 }

  return(obj)

}


plot.rois <- function(obj,which='all',crosshairs=T){

  df <- rois2df(obj)

  if(which !='all')
    df <- subset(df,name %in% which)

  p <- ggplot2::ggplot() + ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)+

    ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),ggplot2::aes(x=xcenter,y=ycenter,label=name),color='white') +
    ggplot2::xlim(c(0,obj$resolution[1])) + ggplot2::ylim(c(obj$resolution[2],0)) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank())


  if(crosshairs){
    p <- p +
    ggplot2::geom_hline(yintercept=384,color='red',size=0.3,linetype='dashed') +
    ggplot2::geom_vline(xintercept=512,color='red',size=0.3,linetype='dashed')
  }

  p
  return(p)

}


 rois2df <- function(obj,which='all'){

  roidfs <- list()

  for(i in 1:length(obj$rois)){

  eltmp <- as.data.frame(vertices(obj$rois[[i]]$roi))
  eltmp$name <- obj$rois[[i]]$name
  eltmp$xcenter <- obj$rois[[i]]$center[1]
  eltmp$ycenter <- obj$rois[[i]]$center[2]
  eltmp$xradius <- obj$rois[[i]]$xradius
  eltmp$yradius <- obj$rois[[i]]$yradius
  eltmp$radius <- obj$rois[[i]]$radius
  eltmp$shape <- obj$rois[[i]]$shape

  roidfs[[i]]<- eltmp

  }

  df <- do.call('rbind',roidfs)
  return(df)


}
