

#class constructor
itrackr <- function(edfs = NULL,path=NULL,pattern=NULL,idvar='ID',resolution=c(1024,768),samples=FALSE)
{

  me <- list(
              edfs = edfs,
              ids = NULL,
              idvar = idvar,
              indexvars = NULL,
              timevars = NULL,
              resolution = resolution,
              header = data.frame,
              samples = data.frame,
              fixations = data.frame,
              saccades = data.frame,
              blinks = data.frame,
              messages = data.frame,
              epochs = list(),
              beh = data.frame
  )


  class(me) <- append(class(me),'itrackR')

  if(!is.null(me$edfs))
    {
      path <- edfs
      pattern <- NULL
    }

  if(!is.null(path))
    me <- load_edfs(me,path=path,pattern=pattern,samples=samples)


  return(me)

}


load_edfs <- function(obj,path='.',pattern='*.edf',samples=FALSE,recursive = FALSE){
  UseMethod('load_edfs',obj)
}

load_edfs.itrackR <- function(obj,path='.',pattern='*.edf',samples=FALSE,recursive = FALSE){

  fields <- c('fixations','saccades','blinks','messages','header')

  if(samples)
    fields <- c(fields,'samples')

  allbatch <- edfR::edf.batch(path,pattern=pattern,samples=samples,do.plot=F)

  alldata <- edfR::combine.eyedata(allbatch,fields=fields)

  obj$header <- alldata$header
  obj$fixations <- alldata$fixations
  obj$saccades <- alldata$saccades
  obj$blinks <-alldata$blinks
  obj$messages <- alldata$messages

  if(samples)
    obj$samples <- alldata$samples

  if(is.null(obj$edfs))
    obj$edfs <- unlist(lapply(tmp,function(x) x$filename))

  obj$edfs <- unlist(lapply(obj$edfs,make_fullpath))


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


plot <- function(obj,zoom=FALSE){
  UseMethod('plot',obj)

}

plot.itrackR <- function(obj,zoom=FALSE){



  plt <- ggplot2::ggplot(obj$fixations,ggplot2::aes(x=gavx,y=gavy)) + ggplot2::geom_point() + ggplot2::facet_wrap(~ID)

  if(zoom)
    plt <- plt+ggplot2::coord_cartesian(xlim=c(0,obj$resolution[1]),ylim=c(0,obj$resolution[2]))

  plt <- plt + ggplot2::theme_bw()

  plt

  return(plt)
}


index_vars <- function(obj,varnames,patterns,numeric.only=FALSE)
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

  obj$header <- dplyr::left_join(obj$header,msg_index,by=c('ID','eyetrial'))

  if(timestamp)
    obj$timevars <- varnames

  return(obj)

}


add_behdata <- function(obj,beh){
  behmerged <- merge(beh,obj$header[c('ID','eyetrial',obj$indexvars)],by=c('ID',obj$indexvars),all.x=T)
  obj$beh <- behmerged
  return(obj)
}


epoch_samples <- function(obj,timevar,field='paL',epoch=c(-100,100))
{

  allids <- unique(obj$header[obj$idvar])

  allepochs <- list()

  for(id in allids)
  {

#
#   thisepoch <- list()
#   thisepoch$name <- timevar
#   thisepoch$epoch <- epoch
#   thisepoch$field <- field

  subdata <- get_subdata(obj,id,fields=c('header','samples','messages'))

  thisepoch <- edfR::epoch.samples(timevar,subdata$samples,sample.field=field,epoch=epoch,eyetrial=T,messages=subdata$messages)

  allepochs <- c(allepochs,thisepoch)
}

  obj$epochs <- allepochs

  return(obj)

}


get_subdata <- function(obj,ID,fields = c('header','samples'))
{
  subobj <- itrackr()

  for(f in fields)
  {

  subobj[[f]] <- obj[[f]][obj[[f]]$ID==ID,]

  }

  return(subobj)




}


load_samples <- function(obj){

  allsamps <- list()

  i <- 1
  for(edf in obj$edfs)
  {
    alldata <- edfR::edf.batch(edf,samples=T,do.plot=F)

    allsamps[[i]] <- alldata$samples
    i <- i+1
  }


  return(allsamps)

}
