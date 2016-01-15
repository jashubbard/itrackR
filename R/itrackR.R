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
              beh = data.frame,
              transform = list()
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
    obj$edfs <- unlist(lapply(allbatch,function(x) x$filename))

    obj$edfs <- unlist(lapply(obj$edfs,make_fullpath))

  return(obj)

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





get_subdata <- function(obj,id,fields = c('header','samples'))
{
  subobj <- itrackr()

  for(f in fields)
  {

    subobj[[f]] <- subset(obj[[f]],ID==id)

  }

  return(subobj)




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

      tmpROI$roi <- spatstat::ellipse(xradius[[i]],yradius[[i]],centre=coords[i,],phi=angles[[i]])


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


eyemerge <- function(obj,eyedata='fixations',behdata='all'){

  eyes <- obj[[eyedata]]
  beh <- obj$beh

  if(behdata !='all')
    beh <- dplyr::select_(beh,.dots=c('ID','eyetrial',behdata))

  output <- dplyr::right_join(beh,eyes,by=c('ID','eyetrial'))
  output <- dplyr::arrange(output,ID,eyetrial)

  return(output)





}


drift_correct <- function(obj,vars=c('ID'),eydata='fixations',threshold = 10){

  fixnames <- names(obj$fixations)

  fixdata <- eyemerge(obj,behdata=unique(c('ID','eyetrial',vars)))

  fixdata <- dplyr::group_by_(fixdata,.dots=vars)
  fixdata <- dplyr::mutate(fixdata,
                              center_x = median(gavx,na.rm=T),
                              center_y = median(gavy,na.rm=T))
  fixdata <- dplyr::ungroup(fixdata)

  fixdata$real_x <- round(obj$resolution[1]/2)
  fixdata$real_y <- round(obj$resolution[2]/2)

  fixdata <- dplyr::mutate(fixdata,
                           gavx_raw = gavx,
                           gavy_raw = gavy,
                           shift_x = center_x - real_x,
                           shift_y = center_y - real_y)


  fixdata$shift_x[abs(fixdata$shift_x)<threshold] <- 0
  fixdata$shift_y[abs(fixdata$shift_y)<threshold] <- 0
  fixdata <- dplyr::mutate(fixdata,
                          gavx = gavx - shift_x,
                          gavy = gavy - shift_y)



 obj$fixations <- fixdata[fixnames]

 obj$transform$fixations <- fixdata[c('ID','eyetrial','fixation_key','shift_x','shift_y')]


return(obj)


}


