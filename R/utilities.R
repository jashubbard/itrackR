
#' Generates a set of "numpoints" coordinates around a circle with specified radius,
#' centered around x and y, and (optionally) with at starting_angle (default=0 degrees)
radialCoords <- function(x,y,numpoints,radius,starting_angle=0){

  coords <- matrix(0,nrow=numpoints,ncol=2)

  for(i in 1:numpoints){

    angle <- (i-1)*(360/numpoints)+starting_angle
    coords[i,1] <- ceiling(x+sin((angle*pi)/180)*radius)
    coords[i,2] <- ceiling(y+cos((angle*pi)/180)*radius)

  }

  colnames(coords) <- c('x','y')

  return(coords)
}

roiFlower <- function(numrois,starting_angle=0,direction='counterclockwise'){

  starting_angle <- pracma::deg2rad(starting_angle)

  angles <- seq(starting_angle,2*pi-starting_angle,(2*pi)/numrois)

  if(direction=='counterclockwise')
    angles <- rev(angles)

  return(angles)
}

calcHits <- function(obj,rois='all',append=FALSE){

  obj <- calcHits_fixations(obj,rois,append)
  obj <- calcHits_saccades(obj,rois,append)

  #record what we did
  funcall <- list(
    rois=rois,
    append=append
  )

  obj <- update_history(obj,'calcHits',funcall,append=F)

  return(obj)
}

calcHits_fixations <- function(obj,rois='all',append=FALSE){

  allrois <- lapply(obj$rois,function(x) x$roi)
  allnames <- unlist(lapply(obj$rois, function(x) x$name))

  allnames <- paste0('roi_',allnames)

  hits <- data.frame(key = obj$fixations$key)

  names(hits)[-1] <- allnames

  for(i in 1:length(allrois)){
    hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$fixations$gavx,obj$fixations$gavy,allrois[[i]]))
  }

  hits <- hits[,c('key',allnames)]

  if(!append)
    tmp <- obj$fixations[c('key','ID','eyetrial','sttime','entime','gavx','gavy')]
  else
    tmp <- obj$fixations

  obj$fixations <- cbind(tmp,dplyr::select(hits,-key))

  return(obj)
}

calcHits_saccades <- function(obj,rois='all',append=FALSE){

  allrois <- lapply(obj$rois,function(x) x$roi)
  allnames <- unlist(lapply(obj$rois, function(x) x$name))

  allnames <- c(paste0('roi_start_',allnames), paste0('roi_end_',allnames))

  hits <- data.frame(key = obj$saccades$key)
  names(hits)[-1] <- allnames

  for(i in 1:length(allrois)){
    hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$saccades$gstx,obj$saccades$gsty,allrois[[i]]))
  }

  for(i in (length(allrois)+1):(length(allrois)*2)){
    hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$saccades$genx,obj$saccades$geny,allrois[[i-length(allrois)]]))
  }

  hits <- hits[,c('key',allnames)]

  if(!append)
    tmp <- obj$saccades[c('key','ID','eyetrial','sttime','entime','gstx','gsty','genx','geny')]
  else
    tmp <- obj$saccades

  obj$saccades <- cbind(tmp,dplyr::select(hits,-key))

  return(obj)
}

mapROIs <- function(obj,names,indicators=NULL,cleanup=FALSE){

  newnames <- names
  startnames <- newnames

  roinames <- unlist(lapply(obj$rois, function(x) x$name))

  for(eyedata in c('fixations','saccades')){

    if(eyedata=='fixations'){
      keyvar='key'
      roivars <- paste0('roi_',roinames)
      newnames <- paste0(startnames,'_hit')
    }
    else if(eyedata=='saccades'){
      keyvar='key'
      roivars = c(paste0('roi_start_',roinames))
      roivars = c(roivars, paste0('roi_end_',roinames))
      newnames <- c(paste0(startnames,'_start_hit'),paste0(startnames,'_end_hit'))
    }

    #if we're re-running this, first remove the existing columns
    if(any(names(obj[[eyedata]]) %in% newnames))
      obj[[eyedata]] <- obj[[eyedata]][ , -which(names(obj[[eyedata]]) %in% newnames)]

    tmp1 <- obj$beh[,c('ID','eyetrial',indicators)] #behavioral data with trial-by-trial indicators
    #tmp1[,indicators] <- as.numeric(tmp1[, indicators]) # if indicators are factors convert them here to numbers
    tmp2 <- obj[[eyedata]][,c('ID','eyetrial',keyvar,roivars)] #fixation data with roi hits
    tmp3 <- dplyr::right_join(tmp1,tmp2,by=c('ID','eyetrial')) #merge them together so we have our indicators for fixation-by-fixation
    tmp3 <- dplyr::arrange_(tmp3,keyvar) #put in the original order
    hits <- tmp3[,roivars] #grab only the hit data (since we will use the indicator to index from this)

    #funny hack to handle situations where indicator is NA.
    #Create a column at the end that is all NAs, so when indexed it fills in that value
    hits$missing <- NA
    hits <- as.matrix(hits)

    for(i in 1:length(indicators)){

      # if indicators are factors / strings (i.e., if ROIs have non-numeric names)
      # map indicator names to corresponding column number

      factorIndicators = c()
      if(!is.numeric(tmp3[,indicators[i]])){
        factorIndicators = tmp3[,indicators[i]]
        tmp3[,indicators[i]] = as.numeric(match(tmp3[,indicators[i]], roinames))
      }

      #wherever indicator is NA, change to #rois+1, so it indexes the "missing" column created above
      tmp3[,indicators[i]][is.na(tmp3[,indicators[i]])] <- ncol(hits)

      tmp3[newnames[i]] <- NA

      idx <- sub2ind(hits,1:nrow(tmp3),tmp3[,indicators[i]]) #get the single number indices based on row and column numbers (Indicator corresponds to column number)
      tmp3[newnames[i]] <- as.numeric(hits[idx]) #grab the actual data based on the single number indices. Yields a single column of hits and misses

      if(eyedata=='saccades'){
        idx2 <- sub2ind(hits,1:nrow(tmp3),tmp3[,indicators[i]]+length(obj$rois)) #get the single number indices based on row and column numbers (Indicator corresponds to column number)
        tmp3[newnames[i+length(startnames)]] <- as.numeric(hits[idx2]) #grab the actual data based on the single number indices. Yields a single column of hits and misses
      }

      tmp3[,indicators[i]][tmp3[,indicators[i]]==ncol(hits)] <- NA #fill the missing data back with NAs
    }

    #restore the factor indicators that were previously replaced with numbers (if appropriate)
    if(length(factorIndicators) != 0){
      tmp3[,indicators[i]] = factorIndicators
    }

    #in case we're re-running this, don't create duplicate columns
    if(any(names(obj[[eyedata]]) %in% indicators)){
      obj[[eyedata]] <- obj[[eyedata]][ , -which(names(obj[[eyedata]]) %in% indicators)]
    }

    #merge our existing fixation data with our new mapped data
    obj[[eyedata]] <- dplyr::left_join(obj[[eyedata]],tmp3[,c(keyvar,indicators,newnames)],by=keyvar)

    # delete the hit columns to free memory
    if(cleanup)
      obj[[eyedata]][, c(roivars)] <- list(NULL)
    # trigger garbage collection:
    gc()


    #keep track of our indicator variables
    if(is.null(obj$roimapvars))
      obj$roimapvars <- indicators
    else
      obj$roimapvars <- unique(c(obj$roimapvars,indicators))

    #record what we did
    funcall <- list()
    funcall$names <- names
    funcall$indicators <- indicators
    funcall$cleanup <- cleanup

    obj <- update_history(obj,'mapROIs',funcall,append=T)



  }
  return(obj)

}

sub2ind <- function(mat, r, c){

  m <- nrow(mat)
  ind <- (c-1)*m + r
  return(ind)
}

ind2sub <- function(mat,ind){
  m <- nrow(mat)
  r = ((ind-1) %% m) + 1
  c = floor((ind-1) / m) + 1
  return(list(row=r,col=c))
}

epoch_fixations <- function(obj,rois=c(1),start=0,end=700,binwidth=25,event='starttime',type='standard',...){

  obj <- epoch_events(obj,eyedata='fixations',rois=rois,epoch = c(start,end),binwidth=binwidth,event=event,type=type,...)
  return(obj)
}

epoch_events <- function(obj,eyedata='fixations',epoch = c(0,700), event='starttime',rois=c(1),binwidth=25,type='standard',saccade_hit = 'end'){

  keyvar <- 'key'
  epochvar <- paste0(eyedata,'_epochs')

  if(eyedata=='fixations')
    saccade_hit <- NULL

  if(eyedata=='blinks')
    rois <- 1

  start <- epoch[1]
  end <- epoch[2]
  startvar <- 'sttime'
  endvar <- 'entime'

  for(roi in rois){

    if(is.numeric(roi))
      hitvar <- paste(c('roi',saccade_hit,roi),collapse = '_')
    else
      hitvar <- paste(c(roi,saccade_hit,'hit'),collapse='_')

    if(eyedata=='blinks'){
      hitvar=NULL
    }

    if(eyedata=='fixations' && !(hitvar %in% names(obj$fixations))){
      warning('running calcHits first...')
      obj <- calcHits(obj,rois=roi,append=T)
    }

    prefix = 't'
    firstbin <- (ceiling(start/binwidth))
    numbins = ((end - start)/binwidth)+1
    lastbin <- firstbin + numbins-1

    if(is.null(event)){
      obj$header$default_event = obj$header$starttime
      event = 'default_event'
    }

    eventdata <- obj$header[,unique(c('ID','eyetrial','starttime',event))]

    df <- obj[[eyedata]][,c('ID','eyetrial',keyvar,startvar,endvar,hitvar)]
    df <- dplyr::left_join(df,eventdata,by=c('ID','eyetrial'))
    df <- dplyr::arrange_(df,.dots = keyvar)

    #convert to "trial time" if we're not time locking to beginning of trial
    if(event !='starttime')
      df[event] <- df[,event] - df$starttime

    df[startvar] <- df[,startvar] - df$starttime
    df[endvar] <- df[,endvar] - df$starttime

    #bin the data based onstart time and bin width (in ms)
    df$bin_start <- ceiling((df[,startvar] - df[,event])/binwidth)

    if(type=='cumulative'){
      df$bin_end <- NA
      df$bin_end[df$bin_start>=firstbin & df$bin_start<=lastbin] <- lastbin
    }
    else
      df$bin_end <- ceiling((df[,endvar] - df[,event])/binwidth)

    #have a column for the bin width (easier for plotting later)
    df$binwidth <- binwidth

    #throw out before stim onset and after numbins
    df <- subset(df,bin_start>=firstbin & bin_end<=lastbin)

    #convert to "trial time"
    #~   df[startvar] <- df[,startvar] - df$starttime
    #~   df[endvar] <- df[,endvar] - df$starttime

    if(start<=0){
      df$start_adj = df$bin_start + abs(firstbin)+1
      df$end_adj = df$bin_end + abs(firstbin)+1
      # intervals <- apply(cbind(start_adj,end_adj),1,function(x) x[1]:x[2])
    }
    else {
      df$start_adj <- df$bin_start
      df$end_adj <- df$bin_end

    }


    # vectors <- intervals2matrix(start_adj,end_adj,df[,hitvar],numbins)

    vars = gsub('-','_',paste0(prefix,seq(firstbin*binwidth,lastbin*binwidth,binwidth)))
    # df[vars] <- vectors

    if(is.null(obj$epochs$fixations[[event]]))
      obj$epochs$fixations[[event]] <- list()

    obj[[epochvar]][[event]][[roi]] <- list()
    obj[[epochvar]][[event]][[roi]]$data <- df[c('ID','eyetrial','key','binwidth','bin_start','bin_end','start_adj','end_adj',hitvar)]
    obj[[epochvar]][[event]][[roi]]$roi <- roi
    obj[[epochvar]][[event]][[roi]]$epoch <- c(start,end)
    obj[[epochvar]][[event]][[roi]]$binwidth <- binwidth
    obj[[epochvar]][[event]][[roi]]$firstbin <- firstbin
    obj[[epochvar]][[event]][[roi]]$lastbin <- lastbin
    obj[[epochvar]][[event]][[roi]]$numbins <- numbins
    obj[[epochvar]][[event]][[roi]]$varnames <- vars

  }

  #record what we did
  funcall <- list()
  funcall$start <- start
  funcall$end <- end
  funcall$binwidth <- binwidth
  funcall$event <- event
  funcall$type <- type
  funcall$roi <- roi

  obj <- update_history(obj,name='epoch_fixations',funcall,append=T)

  return(obj)
}


intervals2matrix <- function(start,end,val,numbins,fillval=NA){


  #new intervals2matrix
  #ge the intervals of start:end bin (1 list per fixation/saccade/whatever)
  intervals <- purrr::map2(start,end,~ seq(.x,.y))
  #figure out the length of each of those intervals
  lengths <- unlist(lapply(intervals,length))
  #we will create a matrix of indices for our big matrix. We need a row num, and for each row, the bin will be the column
  rownums <- unlist(purrr::map2(1:length(lengths),lengths,~ rep(.x,.y)))

  #this creates that matrix. Column 1 is row index, column 2 is column index
  subs <- as.matrix(cbind(rownums,unlist(intervals)))

  #whatever hit variable we want to fill in, get the actual values and replicate them for the full interval
  vals <- unlist(purrr::map2(val, lengths, ~ rep(.x,.y) ))

  #create our sparse matrix, with row and column indices, fill with the values from the hit variable
  sp <- Matrix::sparseMatrix(i=subs[,1],j=subs[,2],x=vals)

  #sparse matrix will have zeros where we didn't specify values
  #for fixations/saccades we want it to be missing data.
  #this creates a "mask" for making those missing values
  if(is.na(fillval)){
    mask <- Matrix::sparseMatrix(i=subs[,1],j=subs[,2],x=rep(999,nrow(subs)))
    nempty <- mask==0
    #cases where fixations/saccades/whatever didn't occur, we make NA
    sp[nempty] <- rep(NA,sum(nempty))
  }else{

    #if it's 0, we don't need this step. Skip to save time
    if(fillval !=0){
      mask <- Matrix::sparseMatrix(i=subs[,1],j=subs[,2],x=rep(999,nrow(subs)))
      nempty <- mask==0
      #cases where fixations/saccades/whatever didn't occur, we make NA
      sp[nempty] <- rep(fillval,sum(nempty))
    }

  }

  return(as.matrix(sp))

}


aggregate_fixation_epochs <- function(obj,event='starttime',rois=c(1),groupvars=c('ID'),level='group',type='probability',condition=NULL,condition.str=FALSE){


  if(!condition.str)
    condition <- deparse(substitute(condition))


  roinames <- list()


  allres <- list()

  for(r in rois){

    df <- eyemerge(obj,'fixations_epochs',event=event,roi=r, behdata = groupvars,condition=condition, condition.str=T)


    timevars <- obj$fixations_epochs[[event]][[r]]$varnames

    if(is.numeric(r))
      hitvar <- paste0('roi_',r)
    else
      hitvar <- paste0(r,'_hit')

    test <- dplyr::select_(df,.dots=paste0('-',hitvar))
    test <- tidyr::gather_(test,'timepoint','val',timevars)
    test <- dplyr::mutate(test,timepoint = gsub('t','',timepoint),
                    timepoint = as.numeric(gsub('_','-',timepoint)))

    if(is.numeric(r)){
      test$roi <- hitvar
      roinames <- c(roinames,hitvar)
    }
    else{
      test$roi <- r
      roinames <- c(roinames,r)
    }

    allres <-c(allres,list(test))
  }


  #aggregate on the trial level
  allres <- dplyr::bind_rows(allres)
  allres <- dplyr::group_by_(allres,.dots=unique(c('ID','eyetrial',groupvars,'timepoint','roi')))
  allres <- dplyr::summarise(allres,val = max(val,na.rm=T))

  if(level %in% c('subject','ID','group')){

    allres <- dplyr::group_by_(allres,.dots=unique(c('ID',groupvars,'timepoint','roi')))
    allres <- dplyr::summarise(allres,val = mean(val,na.rm=T))

    if(tolower(type) %in% c('difference','logratio','proportion')){

      if(length(rois)>2 && type=='difference')
        stop('can only do difference waves with 2 rois at a time')

      allres <- tidyr::spread(allres,roi,val)

      if(tolower(type)=='logratio'){

        allres <- dplyr::mutate_(allres,val = sprintf('log(%s+1) / (%s+1)',roinames[1],roinames[2]),
                                 roi = "'logRatio'")
      }
      else if(tolower(type)=='proportion'){
        allres <- dplyr::mutate_(allres,val = sprintf('%s / %s',roinames[1],roinames[2]),
                                 roi = "'proportion'")
      }
      else if(tolower(type)=='difference'){
        allres <- dplyr::mutate_(allres,val = sprintf('%s - %s',roinames[1],roinames[2]),
                                 roi = "'difference'")
      }

      allres <- dplyr::select_(allres,paste0('-',roinames[1]),
                               paste0('-',roinames[2]))
    }



  }


  if(level=='group'){
    allres <-  dplyr::group_by_(allres,.dots=unique(c(groupvars,'timepoint','roi')))
    allres <- dplyr::summarise(allres, val = mean(val,na.rm=T))
    allres <- dplyr::arrange_(allres,.dots=unique(c(groupvars,'roi','timepoint')))
  }


  return(allres)

}






downsample <- function (v, N){ # v is the input vector, and keep every N sample
  seed <- c(TRUE,rep(FALSE,N-1))
  cont <- rep(seed,ceiling(length(v)/N))[1:length(v)]
  return(v[which(cont)])
}

change.sampledir <- function(obj,dir){

  olddir <- obj$sample.dir

  for(s in obj$samples){

    if(!file.exists(s))
      stop('Some files are missing. Run load.samples to save all sample data first.')

    fname <- basename(s)
    newname <- file.path(dir,fname)
    file.copy(s,newname)
    file.remove(s)
  }

  return(obj)
}

# Implementation of the Engbert & Kliegl algorithm for the
# detection of saccades.  This function takes a data frame of the
# samples and adds three columns:
#
# - A column named "saccade" which contains booleans indicating
#   whether the sample occurred during a saccade or not.
# - Columns named vx and vy which indicate the horizontal and vertical
#   speed.
#code borrowed from the "saccades" package
detect.saccades <- function(samples, lambda=15, smooth.saccades=T) {

  # Calculate horizontal and vertical velocities:
  vx <- stats::filter(samples$x, -1:1/2)
  vy <- stats::filter(samples$y, -1:1/2)

  # We don't want NAs, as they make our life difficult later
  # on.  Therefore, fill in missing values:
  vx[1] <- vx[2]
  vy[1] <- vy[2]
  vx[length(vx)] <- vx[length(vx)-1]
  vy[length(vy)] <- vy[length(vy)-1]

  msdx <- sqrt(median(vx**2, na.rm=T) - median(vx, na.rm=T)**2)
  msdy <- sqrt(median(vy**2, na.rm=T) - median(vy, na.rm=T)**2)

  radiusx <- msdx * lambda
  radiusy <- msdy * lambda

  sacc <- ((vx/radiusx)**2 + (vy/radiusy)**2) > 1
  if (smooth.saccades) {
    sacc <- stats::filter(sacc, rep(1/3, 3))
    sacc <- as.logical(round(sacc))
  }
  samples$saccade <- ifelse(is.na(sacc), F, sacc)
  samples$vx <- vx
  samples$vy <- vy

  samples

}

itrackR.data <- function(data){

  #get path to data folder
  d <- paste0(system.file('extdata', package='itrackR'), .Platform$file.sep)

  #find all edfs, get full path of each
  edfs <-  list.files(path=d,'*.edf')
  edfs <- file.path(d,edfs)

  #return the proper data,depending on what we ask for
  output <- switch(data,
                   path = d,
                   edfs = edfs,
                   beh = readRDS(file.path(d,'beh.rds')))

  return(output)
}


itrackr.data <- function(data){
  output <- itrackR.data(data)
  return(output)
}



remove.subjects <- function(obj,ID){

  if(length(ID)==1)
    ID <- list(ID)

  ID <- ID[ID %in% obj$subs]

  if(length(ID)==0)
    stop("ID(s) not found")

  subs_to_keep <- !(obj$subs %in% ID)

  #remove from our list of edfs and subject ids
  obj$edfs <- obj$edfs[subs_to_keep]
  obj$subs <- obj$subs[subs_to_keep]

  #remove actual data
  obj$header <- obj$header[!(obj$header[[obj$idvar]] %in% ID),]
  obj$fixations <- obj$fixations[!(obj$fixations[[obj$idvar]] %in% ID),]
  obj$saccades <- obj$saccades[!(obj$saccades[[obj$idvar]] %in% ID),]
  obj$blinks <- obj$blinks[!(obj$blinks[[obj$idvar]] %in% ID),]
  obj$messages <- obj$messages[!(obj$messages[[obj$idvar]] %in% ID),]

  #if samples have been loaded
  if(length(z$samples)>0){

    #find the .rds files
    sampfiles <- obj$samples[obj$subs %in% ID]
    obj$samples <- obj$samples[subs_to_keep] #remove from our list

    #delete the actual files
    for(s in sampfiles){
      if(file.exists(s))
        file.remove(s)
    }
  }

  #if drift correction has been performed, remove from the "transform" table
  if(length(obj$transform)>0)
    obj$transform <- obj$transform[!(obj$transform[obj$idvar] %in% ID)]

  #remove from epoched data-- very tedious (rois within events within fixations/samples)

  epoch_types <- names(obj$epochs) #types: fixations/samples

  for(t in epoch_types){
    events <- names(obj$epochs[[t]]) #time-locking events
    for(e in events){
      rois <- names(obj$epochs[[t]][[e]]) #individual rois
      for(r in rois){
        keeprows <- !(obj$epochs[[t]][[e]][[r]][[obj$idvar]] %in% ID) #keep only the rows without the IDs to throw away
        obj$epochs[[t]][[e]][[r]] <- obj$epochs[[t]][[e]][[r]][keeprows,]
      }
    }
  }

  return(obj)
}



add.subjects <- function(obj,edfs,beh){

  comb <- obj

  newsubs <- itrackr(edfs=edfs,resolution=obj$resolution,binocular=obj$binocular)
  newsubs$history <- obj$history

  newsubs <- replay_analysis(newsubs,beh=beh)

  comb$edfs <- c(obj$edfs,newsubs$edfs)
  comb$subs <- c(obj$subs,newsubs$subs)
  comb$fixations <- rbind(obj$fixations,newsubs$fixations)
  comb$fixations$key <- 1:nrow(comb$fixations)

  comb$saccades <- rbind(obj$saccades,newsubs$saccades)
  comb$saccades$key <- 1:nrow(comb$saccades)

  comb$messages <- rbind(obj$messages,newsubs$messages)
  comb$messages$key <- 1:nrow(comb$messages)

  comb$blinks <- rbind(obj$blinks,newsubs$blinks)
  comb$blinks$key <- 1:nrow(comb$blink)

  comb$beh <- rbind(obj$beh,newsubs$beh)

  #if drift correction has been performed, remove from the "transform" table
  if(length(obj$transform)>0)
    comb$transform <- rbind(obj$transform,newsubs$transform)

  #remove from epoched data-- very tedious (rois within events within fixations/samples)

  epoch_types <- names(obj$epochs) #types: fixations/samples

  for(t in epoch_types){
    events <- names(obj$epochs[[t]]) #time-locking events
    for(e in events){
      rois <- names(obj$epochs[[t]][[e]]) #individual rois
      for(r in rois){

        comb$epochs[[t]][[e]][[r]] <- rbind(obj$epochs[[t]][[e]][[r]],newsubs$epochs[[t]][[e]][[r]])
      }
    }
  }


  return(comb)
}



reset <- function(obj, reload = FALSE, rebuild = FALSE,keep.rois = FALSE, beh=NULL){


  if(reload){
    newobj <- itrackr(edfs = obj$edfs, resolution = obj$resolution, binocular = obj$binocular)

    if(keep.rois)
      newobj$rois <- obj$rois
  }
  else{

    #undo drift correction if it was done
    obj <- undrift(obj)

    newobj <- itrackr() #blank itrackr object


    newobj$fixations <- dplyr::select(obj$fixations, -dplyr::matches('roi_*'),-dplyr::matches('*_hit'))
    newobj$fixations <- newobj$fixations[!(names(newobj$fixations) %in% obj$roimapvars)]
    newobj$saccades <- dplyr::select(obj$saccades, -dplyr::matches('roi_*'),-dplyr::matches('*_hit'))
    newobj$saccades <- newobj$saccades[!(names(newobj$saccades) %in% obj$roimapvars)]
    newobj$header <- obj$header[!(names(obj$header) %in% c(obj$roimapvars,obj$indexvars, obj$timevars))]
    newobj$edfs <- obj$edfs
    newobj$ids <- obj$ids
    newobj$binocular <- obj$binocular
    newobj$resolution <- obj$resolution
    newobj$sample.dir <- obj$sample.dir
    newobj$idvar <- obj$idvar
    newobj$blinks <- obj$blinks
    newobj$messages <- obj$messages
    newobj$samples <- obj$samples

    if(keep.rois)
      newobj$rois <- obj$rois

  }

  if(rebuild){
    print('rebuilding itrackr object using history of analysis')
    newobj$rois <- list()
    newobj$history <- obj$history

    if(is.null(beh))
      newobj$beh <- obj$beh[!(names(obj$beh) %in% c(obj$timevars))]

    newobj <- replay_analysis(newobj,beh=beh)
  }

  return(newobj)
}

