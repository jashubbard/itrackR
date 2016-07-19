#class constructor
itrackr <- function(txt = NULL,edfs = NULL,path=NULL,pattern=NULL,resolution=c(1024,768),datadir=tempdir(),binocular=FALSE)
{

  me <- list(
              edfs = edfs,
              ids = NULL,
              idvar = 'ID',
              binocular = binocular,
              indexvars = NULL,
              timevars = NULL,
              resolution = resolution,
              header = data.frame(),
              samples = list(),
              sample.dir = NULL,
              fixations = data.frame(),
              saccades = data.frame(),
              blinks = data.frame(),
              messages = data.frame(),
              epochs = list(),
              beh = data.frame(),
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

  if(!is.null(txt))
    me <- load_txt(me,filename=txt)

  return(me)

}

load_edfs <- function(obj,path='.',pattern='*.edf',recursive = FALSE){

  fields <- c('fixations','saccades','blinks','messages','header')


  allbatch <- edfR::edf.batch(path,pattern=pattern,samples=F,do.plot=F)

  alldata <- edfR::combine.eyedata(allbatch,fields=fields)

 #create key variables for each, for easy sorting
  alldata$fixations$fixation_key <- 1:nrow(alldata$fixations)
  alldata$saccades$saccade_key <- 1:nrow(alldata$saccades)
  alldata$blinks$blink_key <- 1:nrow(alldata$blinks)
  alldata$messages$message_key <- 1:nrow(alldata$messages)


  #get the baseline for each subject, repeated for all rows of their data in the header
  header <- dplyr::group_by_(alldata$header, .dots=(obj$idvar)) %>%
    dplyr::mutate(first_sample = first(starttime)-1) %>% dplyr::ungroup(.)
  obj$header <- timeshift(header, header$first_sample, c('starttime','endtime')) #shift by those times

  #now get a table of IDs and first_sample for shifting everything else
  firsts <- dplyr::distinct_(header,.dots=c(obj$idvar,'first_sample'))

  #loop through all the tables
  for(tbl in c('fixations','saccades','blinks','messages')){

    tmp <- alldata[[tbl]]
    tmp$rownum <- 1:nrow(tmp) #this is so we can re-order it easily

    tmp <- dplyr::left_join(tmp,firsts,by=obj$idvar) #merge with the firsts table
    tmp <- dplyr::arrange(tmp,rownum) #order it back to original
    tmp <- dplyr::select(tmp,-rownum) #throw away our order variable

    if(tbl=='messages')
      vars_to_shift = 'sttime'
    else
      vars_to_shift = c('sttime','entime')

    #time-shift the appropriate columns
    obj[[tbl]] <- timeshift(alldata[[tbl]],tmp$first_sample,vars_to_shift)

  }



#   if(samples)
#     {
#       allsamples <- load_samples(obj)
#       obj$samples <- allsamples$samples
#     }

  if(is.null(obj$edfs))
    obj$edfs <- unlist(lapply(allbatch,function(x) x$filename))

    obj$edfs <- unlist(lapply(obj$edfs,make_fullpath))


  obj$subs <- edf2id(obj$edfs)
  return(obj)

}


load_txt <- function(obj,type='fixations',filename,sep='\t'){


if(type=='fixations')
  fnames <- c('RECORDING_SESSION_LABEL','TRIAL_START_TIME','CURRENT_FIX_START','CURRENT_FIX_END','CURRENT_FIX_X','CURRENT_FIX_Y','TRIAL_INDEX')
if(type=='saccades')
  fnames <- c('RECORDING_SESSION_LABEL','TRIAL_START_TIME','CURRENT_SAC_START_TIME','CURRENT_SAC_END_TIME','CURRENT_SAC_START_X','CURRENT_SAC_START_Y',
              'CURRENT_SAC_END_X','CURRENT_SAC_END_Y')
if(type=='messages')
  fnames = c()

fixdata <- read.delim(filename,sep=sep)

if(any(!(fnames %in% names(fixdata)))){
  missing_names <- fnames[!fnames %in% names(fixdata)]
  error(paste0("The following fields are needed in the input file: ",toString(fnames)))
}

if(type=='fixations'){
fixdata <- dplyr::mutate(fixdata,
                         ID = RECORDING_SESSION_LABEL,
                         eyetrial = TRIAL_INDEX,
                         sttime = TRIAL_START_TIME + CURRENT_FIX_START,
                         entime = TRIAL_START_TIME + CURRENT_FIX_END,
                         gavx = CURRENT_FIX_X,
                         gavy = CURRENT_FIX_Y,
                         fixation_key = 1:nrow(fixdata))
}
else if(type=='saccades'){


}

fixdata$ID <- as.numeric(gsub("([0-9]*).*","\\1",fixdata$ID))

obj$fixations <- dplyr::select(fixdata,ID,eyetrial,sttime,entime,gavx,gavy,fixation_key)

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
#~     eyedata <- cbind(eyedata,obj$header[obj$timevars] - obj$header$starttime)
    eyedata <- cbind(eyedata,obj$header['starttime'])

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


# coords defaults to (0, 0) so you can specify a polygon without the need to specify coords (it is not used because all coordinates are inside the poylgon data frame)
makeROIs <- function(obj,coords=data.frame(x=c(0),y=c(0)),shapes='circle',radius=0,xradius=0,yradius=0,angles=NULL,names=NULL,polygon=data.frame(x=c(),y=c()),append=F){

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

  if(is.null(names)){
    names <- start_roi:(start_roi-1+nrow(coords))
  }

	roipos <- start_roi

	if(shapes[[1]] == 'polygon'){
		# for a polygon only one ROI at a time is possible
			tmpROI <- list()

			tmpROI$roi = spatstat::owin(poly=polygon)
			tmpROI$name <- names[[1]]
			tmpROI$shape <- 'polygon'
			tmpROI$center <- spatstat::centroid.owin(tmpROI$roi)
			tmpROI$radius <- NA
			tmpROI$xradius <- NA
			tmpROI$yradius <- NA

			obj$rois[[roipos]] = tmpROI
	} else {

		for(i in 1:nrow(coords)){

			tmpROI <- list()

			if(shapes[[i]]=='circle'){
				tmpROI$roi <- spatstat::disc(radius=radius[[i]],centre=coords[i,])
			}
			else if(shapes[[i]]=='ellipse'){
				tmpROI$roi <- spatstat::ellipse(xradius[[i]],yradius[[i]],centre=coords[i,],phi=angles[[i]])
			}

			tmpROI$name <- names[[i]]
			tmpROI$shape <- shapes[[i]]
			tmpROI$center <- coords[i,]
			tmpROI$radius <- radius[[i]]
			tmpROI$xradius <- xradius[[i]]
			tmpROI$yradius <- yradius[[i]]

			obj$rois[[roipos]] <- tmpROI
			roipos <- roipos + 1
		}
	}
  return(obj)

}


eyemerge <- function(obj,eyedata='fixations',behdata='all',all.rois=F,event=NULL,roi=NULL,trialtime = TRUE,condition=NULL,condition.str=FALSE){

  #if we're getting epoched fixation data, we already have what we need, just need to grab the right variables, etc.
  if(eyedata=='epoched_fixations'){

    if(is.null(event) || is.null(roi))
      stop('You must provide the name of the time-locking event and the ROI')

    eyes <- obj$epochs$fixations[[event]][[roi]]
    beh <- obj$beh

    realbehvars <- setdiff(colnames(beh),c('ID','eyetrial',obj$indexvars))

    if(behdata[1] !='all')
      behvars <- intersect(behdata,realbehvars)
    else
      behvars <- intersect(colnames(beh),realbehvars)

    timevars <- names(eyes)[grepl('^t[1-9]',names(eyes))]

    eyevars <- c('ID','eyetrial',obj$indexvars,'roi','epoch_start','epoch_end','binwidth','sttime','entime',paste0(roi,'_hit'),timevars)
    eyes <- eyes[eyevars]

    output <- dplyr::right_join(beh,eyes,by=c('ID','eyetrial',obj$indexvars))
    output <- dplyr::arrange(output,ID,eyetrial)
  }
  else{

    #regular fixations or saccades
    eyes <- obj[[eyedata]]
    beh <- obj$beh

    hdr <- obj$header

    #we want some info from the header. Only any events that we found.
    event_vars <- names(hdr)[-which(names(hdr) %in% c('ID','eyetrial','starttime','endtime','duration','first_sample', obj$indexvars))]

    if(length(event_vars)>0){
      hdr <- hdr[c('ID','eyetrial',obj$indexvars,event_vars)]
      beh <- dplyr::left_join(beh,hdr,by=c('ID','eyetrial'))
    }
    #add trial header information
    # beh <- dplyr::left_join(beh,obj$header[c('ID','eyetrial','starttime')], by=c('ID','eyetrial'))


    #remove the roi_1, roi_2, ... columns
    if(!all.rois && any(grepl('^roi_*',colnames(eyes)))){
      eyes <- dplyr::select(eyes,-matches("^roi_*"))
    }

    # #if we want only some of the behavioral variables
    # if(behdata[1] !='all')
    #   beh <- dplyr::select_(beh,.dots=unique(c('ID','eyetrial','starttime',obj$indexvars,behdata)))

    if(behdata[1]=='all')
      behvars <- names(beh)
    else{
      realbehvars <- setdiff(colnames(beh),c('ID','eyetrial',obj$indexvars))
      behvars <- intersect(colnames(beh),realbehvars)
    }

    #in case there are variable names in common between eyedata and behdata, besides index variables
    #remove them from eyedata
    realvars <- setdiff(colnames(eyes),c('ID','eyetrial','starttime',obj$indexvars))
    commonvars <- intersect(realvars,colnames(beh))

    if(length(commonvars)>0)
      eyevars <- names(eyes)[!(colnames(eyes) %in% commonvars)]
      eyes <- eyes[eyevars]

    #merge behavioral and eye data
    output <- dplyr::right_join(beh,eyes,by=c('ID','eyetrial'))
    output <- dplyr::arrange(output,ID,eyetrial)

    output <- dplyr::rename(output,
                           trialsttime = starttime)

    #by default, change starting and endind times of fixations relative to trial start
    if(trialtime){
      output$sttime <- output$sttime - output$trialsttime
      output$entime <- output$entime - output$trialsttime

      if(length(event_vars)>0)
        output[event_vars] <- output[event_vars] - output$trialsttime

    }
  }

  #filter out data based on some condition
  if(!condition.str)
    condition <- deparse(substitute(condition))

   if(condition!="NULL"){
    condition_call <- parse(text=condition)
    r <- eval(condition_call,output, parent.frame())
    output <- output[r, ]
}

  #keep only some behavioral variables
  if(behdata[1]!='all')
    output <- dplyr::select_(output, .dots=unique(c('ID','eyetrial',behvars,names(eyes))))


  return(output)
}

drift_correct <- function(obj,vars=c('ID'),eydata='fixations',threshold = 10){

  fixnames <- names(obj$fixations)

  fixdata <- eyemerge(obj,behdata=unique(c('ID','eyetrial',vars)),trialtime=FALSE)

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
