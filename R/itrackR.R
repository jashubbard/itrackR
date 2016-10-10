#' @title Create an itrackR object
#'
#' @description
#' Create an itrackR object by loading edfs
#'
#' @param edfs list of edf file names as strings. Not necessary if you use \code{path} and \code{pattern}
#' @param path if not providing edf file names, you can provide a path to a directory and a search string
#' @param pattern search string for finding edf files, as a regular expression (default = '*.edf')
#' @param resolution a list specifying the x and y resolution of the monitor (default = c(1024,768))
#' @param binocular whether the experiment used binocular tracking (currently not supported, default = FALSE)

#' @author Jason Hubbard, \email{hubbard3@@uoregon.edu}
#'
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Loads all edf files starting with 3 from a certain directory
#' z <- itrackr(path='/path/to/some/directory', pattern = '3*.edf')
#' }
#'
itrackr <- function(edfs = NULL,path=NULL,pattern='*.edf',resolution=c(1024,768),binocular=FALSE)
{

  me <- list(
              edfs = edfs,
              ids = NULL,
              idvar = 'ID',
              binocular = binocular,
              indexvars = NULL,
              roimapvars = NULL,
              timevars = NULL,
              resolution = resolution,
              header = data.frame(),
              samples = list(),
              sample.dir = NULL,
              fixations = data.frame(),
              saccades = data.frame(),
              blinks = data.frame(),
              messages = data.frame(),
              rois = list(),
              epochs = list(),
              beh = data.frame(),
              transform = list(),
              history = list(step=0)
  )


  class(me) <- append('itrackR',class(me))

  if(!is.null(me$edfs))
    {
      path <- edfs
      pattern <- NULL
    }

  if(!is.null(path))
    me <- load_edfs(me,path=path,pattern=pattern)

  # if(!is.null(txt))
  #   me <- load_txt(me,filename=txt)

  return(me)

}

#' @title load edf files
#'
#' @description
#' Internal function used by \code{itrackr} to load edfs. Do not call directly.
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


  if(is.null(obj$edfs))
    obj$edfs <- unlist(lapply(allbatch,function(x) x$filename))

    obj$edfs <- unlist(lapply(obj$edfs,make_fullpath))


  obj$subs <- edf2id(obj$edfs)
  return(obj)

}

#' @title load from text file from Eyelink DataViewer reports
#'
#' @description
#' Experimental function for loading data from text files created from reports in DataViewer. Currently unused.
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


#' @title Find messages to use as index variables in itrackR object
#'
#' @description
#' Mostly a wrapper for \link{\code{find_messages}}, but extracts messages and sets them as \code{obj$indexvars}. These are used in other functions
#' like \link{\code{add_behdata}} and \link{\code{eyemerge}} to use as index variables for merging with behavioral data. This assumes you sent messages
#' to Eyelink during an experiment to specify the trial number. It can also handle multiple variables (e.g., Block 5, Trial 6).
#'
#' @param obj an itrackR object
#' @param varnames a list specifying what we want the variable names to be called
#' @param patterns search strings for the messages of interest (can use regular expressions)
#' @param numeric.only return only numeric information from the messages (default = FALSE). Useful for indexing variables (e.g., "Trial 5" is converted to 5).
#'
#'
#' @return all variables are stored in \code{obj$header}, based on the \code{varnames} specified. You can also see them in \code{obj$indexvars}.
#' When calling \code{eyemerge}, these will be appended to the data frame. When using \code{\link{add_behdata}}, variable names
#' in the behavioral file should match \code{varnames}.
#'
#'
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
#' # 1 and 7.
#' z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*'), numeric.only = TRUE)
#'
#' #you can see the result by looking at the header
#' View(z$header)
#' }
#'
#'
#' @seealso \code{\link{add_behdata}}, \code{\link{find_messages}}
#'
set_index <- function(obj,varnames,patterns,numeric.only=FALSE)
{
  obj <- find_messages(obj,varnames,patterns,numeric.only)

  obj$indexvars <- varnames

  #record what we did for later

  funcall <- list()
  funcall$varnames <- varnames
  funcall$patterns <- patterns
  funcall$numeric.only = numeric.only

  obj <- update_history(obj,'set_index',funcall,append=F)


  return(obj)
}

#' @title Find messages sent to Eyelink
#'
#' @description
#' General function for finding messages sent to Eyelink and extracting information from them. Useful if you want timestamps
#' of particular events merged with behavioral data. Also used by \code{\link{set_index}} to extract information for index variables (for merging).
#'
#' @param obj an itrackR object
#' @param varnames a list specifying what we want the variable names to be called
#' @param patterns search strings for the messages of interest (can use regular expressions)
#' @param numeric.only return only numeric information from the messages (default = FALSE). Useful for indexing variables (e.g., "Trial 5" is converted to 5).
#' @param timestamp return the timestamp of the message, not the message itself (default = FAlSE)
#'
#'
#' @return all messages are stored in \code{obj$header}, based on the \code{varnames} specified.
#' When calling \code{eyemerge}, these will be appended to the data frame.
#'
#'
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Loads all edf
#' z <- find_messages(z,c("STIMONSET"), patterns=c('STIMONSET'), timestamp = TRUE)
#' }
#'
#'
#' @seealso \code{\link{set_index}}
#'
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

  #record what we did
  obj$history$step <- obj$history$step+1

  funcall <- list()
  funcall$varnames <- varnames
  funcall$patterns <- patterns
  funcall$numeric.only <- numeric.only
  funcall$timestamp <- timestamp

  obj <- update_history(obj,'find_messages',funcall,append=T)

  if(timestamp)
    obj$timevars <- varnames

  return(obj)

}


#' @title Merge behavioral data with itrackR object
#'
#' @description
#' Take a data frame containing behaivoral data and merge with the eyetracking data in an itrackR object. Merging is done based on \code{ID} and any
#' index variables set using \code{\link{set_index}}. Index variables should match the same ones that are in \code{obj$header}.
#'
#'
#'
#' @param obj an itrackR object
#' @param beh a data frame containing behavioral data. Should at least have the column \code{ID} which matches ID names/numbers in \code{obj$subs}.
#' @param append set to TRUE if you've already added behavioral data and you're adding more columns (default = FALSE, overwrite existing data).

#'
#'
#' @return data frame is saved in \code{obj$beh}, with the extra column \code{eyetrial} which matches the eyetracking data.
#' Variables in \code{obj$beh} can be used for subsetting when calling \code{\link{eyemerge}}, \code{\link{plot.itrackR}}, and \code{\link{fixation_timeseries}}.
#'
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
#' # 1 and 7.
#' z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*), numeric.only = TRUE)
#'
#' beh <- itrackr.data('beh)
#'
#' z <- add_behdata(z,beh)
#'
#' #notice the new eyetrial variable...
#' View(z$beh)
#'
#' }
#'
#'
#' @seealso \code{\link{eyemerge}}
#'
add_behdata <- function(obj,beh=NULL,append=FALSE){

  if(append){

    obj$beh <- merge(obj$beh,beh,by=c('ID',obj$indexvars))

    }
  else{

    eyedata <- obj$header[c('ID','eyetrial',obj$indexvars,'starttime')]
    #~     eyedata <- cbind(eyedata,obj$header[obj$timevars] - obj$header$starttime)
    # eyedata <- cbind(eyedata,obj$header['starttime'])

    if(!is.null(beh)){
      behmerged <- merge(beh,eyedata,by=c('ID',obj$indexvars),all.x=T)
      obj$beh <- behmerged
    }
    else
      obj$beh <- eyedata

    }

  #record what we did
  funcall <- list()
  funcall$behnames <- names(beh)

  obj <- update_history(obj,'add_behdata',funcall,append=F)

  return(obj)
}

#' @title Create regions of interest (ROIs) for itrackR object
#'
#' @description
#' Create regions of interest (ROIs) for an itrackR object. ROIs can be circular, elliptical, or polygons. After creating ROIs, you use \code{\link{calcHits}}
#' to determine whether each fixation/saccade fell within each ROI. The helper functions \code{\link{radialCoords}} and \code{\link{roiFlower}} can be used to
#' generate coordinates around some set of coordinates.
#'
#'
#' @param obj an itrackR object
#' @param coords data frame or matrix specifying x and y coordinates for the center of each ROI (default = (0,0), upper-left corner)
#' @param shapes shape for each ROI listed, can be \code{'circle'}, \code{'ellipse'}, or \code{'polygon'}
#' @param radius for circular ROIs, radius of the circle (in pixels)
#' @param xradius for elliptical ROIs, radius along the x dimension
#' @param yradius for elliptical ROIs, radius along the y dimension
#' @param angles for elliptical ROIs, rotation for each ROI in degrees.
#' @param names optional, name for each ROI. Defaults to 1 to number of ROIs listed. ROIs should not have duplicate names!
#' @param polygon for polygonal ROIS, data frame or matrix specifying corners of polygon
#' @param append whether we're adding to our list of ROIs or overwriting (default = FALSE, overwrite)

#'
#'
#' @return ROI information is stored in \code{obj$rois}. Uses the \code{\link{spatstat}} package to store ROI information.
#' To get ROI information as a data frame, you can use the internal function \code{\link{roi2df}}. ROIs can be plotted using \code{\link{plot.rois}}.
#'
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' coords <- matrix(c(512,384,100,100),nrow=2,byrow=T)
#'
#' z <- makeROIs(z,coords,shape='ellipse',xradius=60, yradius=120, angles = c(45, 75))
#' }
#'
#'
#' @seealso  \code{\link{plot.rois}} \code{\link{radialCoords}}  \code{\link{roiFlower}} \code{\link{calcHits}}
#'
makeROIs <- function(obj,coords=NULL,shapes='circle',radius=100,xradius=100,yradius=150,angles=NULL,
                     names=NULL,polygon=data.frame(x=c(),y=c()),append=F){

  if(is.null(coords)){
    #default to center of the screen
    coords <- matrix(round(obj$resolution/2),nrow=1)

  }

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

	#record what we did
	funcall <- list(coords=coords,
	                shapes=shapes,
	                radius=radius,
	                xradius=xradius,
	                yradius=yradius,
	                angles=angles,
	                names=names,
	                polygon=polygon,
	                append=append)

	obj <- update_history(obj,'makeROIs',funcall,append=T)


  return(obj)

}


#' @title Create data frame of merged eyetracking and behavioral data
#'
#' @description
#' Function for creating a flat table of eyetracking and behavioral data from itrackR object. Can be used for gathering fixation, saccade or blink data,
#' as well as fixation data epoched around some event. Will include all behavioral variables, or a subset of named variables.
#'
#'
#'
#' @param obj an itrackR object
#' @param eyedata type of eyetracking data to include. Can be:
#' \itemize{
#'  \item \code{'fixations'} (default). Fixation start/end times and x/y coordinates
#'  \item \code{'saccades'} Saccade start/end times and star/end x/y coordinates
#'  \item \code{'blinks'} Blink start/end times
#'  \item \code{'epoched_fixations'} fixation data epoched around some event using \code{\link{epoch_fixations}}
#' }
#'
#' @param behdata list of variables from \code{obj$beh} to include in the output (default = 'all').
#' @param all.rois whether to include hits/misses for all rois in itrackR object (default = FALSE).
#' @param event the timelocking event when getting epoched fixations
#' @param roi the roi when getting epoched fixations
#' @param trialtime whether to convert fixation/saccade times relative to the current trial, or leave relative to whole experiment.
#' @param condition subset data according to some behavioral variables. Can use unquoted variable names as in \code{\link{subset}}. (e.g., condition = Block <= 5)
#' @param condition.str whether the condition is saved as a string ('Block <= 5'). Used internally.
#'
#' @return returns a data frame with the variables from \code{obj$beh} merged with the eyetracking data requested. Merging is based on \code{ID}, and \code{obj$indexvars}.
#' This corresponds to \code{obj$beh$eyetrial} as well.
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
#' # 1 and 7.
#' z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*), numeric.only = TRUE)
#'
#' beh <- itrackr.data('beh)
#'
#' z <- add_behdata(z,beh)
#'
#' #get fixation data
#' fixes <- eyemerge(z,'fixations')
#'
#' #get data only when Task==2
#' fixes <- eyemerge(z, 'fixations', condition = Task==2)
#'
#' #saccades
#' saccs <- eyemerge(z,'saccades')
#'
#' }
#'
#'
#' @seealso \code{\link{epoch_fixations}} \code{\link{set_index}}
#'
eyemerge <- function(obj,eyedata='fixations',behdata='all',all.rois=F,event='starttime',roi=NULL,trialtime = TRUE,condition=NULL,condition.str=FALSE){

  if(length(obj$beh)==0){
    warning('No behavioral data has been added. Creating dummy behavioral data')
    obj <- add_behdata(obj)
  }

  #if we're getting epoched fixation data, we already have what we need, just need to grab the right variables, etc.
  if(eyedata=='epoched_fixations'){

    if(is.null(event) || is.null(roi))
      stop('You must provide the name of the time-locking event and the ROI')

    eyes <- obj$epochs$fixations[[event]][[roi]]
    beh <- obj$beh

    realbehvars <- setdiff(colnames(beh),c('ID','eyetrial',obj$indexvars,event))

    if(behdata[1] !='all')
      behvars <- intersect(behdata,realbehvars)
    else
      behvars <- intersect(colnames(beh),realbehvars)

    timevars <- names(eyes)[grepl('^t[_1-9]',names(eyes))]

    if(is.numeric(roi))
      hitvar <- paste0('roi_',roi)
    else
      hitvar <- paste0(roi,'_hit')


    eyevars <- c('ID','eyetrial',obj$indexvars,'roi','epoch_start','epoch_end','binwidth','sttime','entime',hitvar,timevars)
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
      hdr <- hdr[c('ID','eyetrial',event_vars)]
      beh <- dplyr::left_join(beh,hdr,by=c('ID','eyetrial'))
    }
    #add trial header information
    # beh <- dplyr::left_join(beh,obj$header[c('ID','eyetrial','starttime')], by=c('ID','eyetrial'))


    #remove the roi_1, roi_2, ... columns
    if(!all.rois && any(grepl('^roi_*',colnames(eyes)))){
      eyes <- dplyr::select(eyes,-dplyr::matches("^roi_*"))
    }

    # #if we want only some of the behavioral variables
    # if(behdata[1] !='all')
    #   beh <- dplyr::select_(beh,.dots=unique(c('ID','eyetrial','starttime',obj$indexvars,behdata)))

    if(behdata[1]=='all')
      behvars <- names(beh)
    else{
      # realbehvars <- setdiff(colnames(beh),c('ID','eyetrial',obj$indexvars))
      # behvars <- intersect(colnames(beh),realbehvars)

      behvars <- unique(c('ID','eyetrial',obj$indexvars,behdata))
    }

    #in case there are variable names in common between eyedata and behdata, besides index variables
    #remove them from eyedata
    realvars <- setdiff(colnames(eyes),c('ID','eyetrial','starttime',obj$indexvars))
    commonvars <- intersect(realvars,colnames(beh))

    if(length(commonvars)>0){
      eyevars <- names(eyes)[!(colnames(eyes) %in% commonvars)]
      eyes <- eyes[eyevars]
    }
    #merge behavioral and eye data
    output <- dplyr::right_join(beh,eyes,by=c('ID','eyetrial'))
    output <- dplyr::arrange(output,ID,eyetrial)



    #by default, change starting and endind times of fixations relative to trial start
    if(trialtime){
      output$sttime <- output$sttime - output$starttime
      output$entime <- output$entime - output$starttime

      if(length(event_vars)>0)
        output[event_vars] <- output[event_vars] - output$starttime

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
    output <- dplyr::select_(output, .dots=unique(c('ID','eyetrial',obj$indexvars,behvars,'starttime',names(eyes))))

  #to make things less confusing
  output <- dplyr::rename(output,
                          trialsttime = starttime)

  return(output)
}



#' @title Perform drift correction on fixations for itrackR object
#'
#' @description
#' Correct for drifts in fixation data. Algorithm computes the median x/y coordinates separately for each subject and behavioral variable specified,
#' then fixations are adjusted based on difference between this "true" center and the center of the screen based on \code{obj$resolution}. Fixations are only
#' adjusted if deviation is greater than some threshold.
#'
#'
#' @param obj an itrackR object
#' @param vars variables to specify subsets of data for performing drift correction. default is \code{'ID'}. Can also name any column in \code{obj$beh} (e.g., 'Block')
#' @param append set to TRUE if you've already added behavioral data and you're adding more columns (default = FALSE, overwrite existing data).
#' @param threshold do not adjust fixations that deviate less than this threshold (in pixels). Default is 10 pixels
#'
#'
#' @return All fixations in \code{obj$fixations} are adjusted, and adjustment amounts are stored in \code{obj$transform}.
#' To undo correction, use \code{\link{undrift}}.
#
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
#' # 1 and 7.
#' z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*), numeric.only = TRUE)
#'
#' beh <- itrackr.data('beh)
#'
#' z <- add_behdata(z,beh)
#'
#' #perform drift correction for each subject and Block.
#' z <- drift_correct(z, vars=c('Block'))
#'
#' #if you're curious how adjustment was done, check it out:
#'
#' View(z$transform)
#'
#' }
#'
#'
#' @seealso \code{\link{undrift}}
#'
drift_correct <- function(obj,vars=c('ID'),threshold = 10){

  fixnames <- names(obj$fixations)

  fixdata <- eyemerge(obj,behdata=unique(c('ID','eyetrial',vars)),trialtime=FALSE) %>%
    dplyr::arrange(fixation_key)

  fixdata <- dplyr::group_by_(fixdata,.dots=unique(c('ID',vars))) %>%
    dplyr::mutate(center_x = median(gavx,na.rm=T),
                  center_y = median(gavy,na.rm=T)) %>%
    dplyr::ungroup(.)


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
                          gavy = gavy - shift_y) %>%
    dplyr::arrange(fixation_key)


  #add the new gavx and gavy values
  #merge so it's robust to weird quirks when calling eyemerge (duplicate rows, etc)
  obj$fixations <- dplyr::left_join(dplyr::select(obj$fixations,-gavx,-gavy),dplyr::select(fixdata,fixation_key,gavx,gavy),by='fixation_key') %>%
    dplyr::arrange(fixation_key)


  # obj$fixations$gavx <- fixdata$gavx
  # obj$fixations$gavy <- fixdata$gavy

  obj$transform$fixations <- fixdata[c('ID','eyetrial','fixation_key','shift_x','shift_y')]

  #record what we did

  funcall <- list()
  funcall$eyedata <- fixdata
  funcall$vars <- vars
  funcall$threshold <- threshold

  obj <- update_history(obj,'drift_correct',funcall,append=F)


return(obj)
}


#' @title Reverse drift correction performed on an itrackR object
#'
#' @description
#' In case you want to undo the drift correction performed by \code{\link{drift_correct}}. It puts fixations back to original coordinates and deletes \code{obj$transform}.
#'
#'
#' @param obj an itrackR object
#'
#' @return All fixations in \code{obj$fixations} are adjusted back to original coordinates. \code{obj$transform} is deleted.

#
#' @examples
#' \dontrun{
#' # itrackr.data('edfs') returns full path to 2 edf files
#' z <- itrackr(edfs=itrackr.data('edfs'))
#'
#' #Extracts messages like "BLOCK 1", and "TRIAL 7" and creates variables "Block" and "Trial", giving them values
#' # 1 and 7.
#' z <- set_index(z,c('Block', 'Trial'), patterns=c('BLOCK [0-9]*', 'TRIAL [0-9][0-9]*), numeric.only = TRUE)
#'
#' beh <- itrackr.data('beh)
#'
#' z <- add_behdata(z,beh)
#'
#' #fixations before correcting
#' plot(z)
#'
#' #perform drift correction for each subject and Block.
#' z <- drift_correct(z, vars=c('Block'))
#'
#' #after correcting
#' plot(z)
#'
#'
#' #never mind
#' z <- undrift(z)
#'
#' #back to normal
#' plot(z)
#'
#' }
#'
#'
#' @seealso \code{\link{drift_correct}}
#'
undrift <- function(obj){

  #if there's no drift correction, just go back
  if(length(obj$transform)==0)
    return(obj)

  fixdata <- dplyr::left_join(obj$fixations, obj$transform$fixations,by=c(obj$idvar,'eyetrial','fixation_key'))

  fixdata <- dplyr::mutate(fixdata,
                           gavx = gavx + shift_x,
                           gavy = gavy + shift_y) %>%
    dplyr::arrange(fixation_key)

  obj$fixations <- fixdata[names(obj$fixations)]

  obj$transform$fixations <- NULL

  #erase the record!
  obj$history$step <- obj$history$step -1
  obj$history$drift_correct <- NULL

  return(obj)

}
