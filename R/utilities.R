
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

  return(obj)
}

calcHits_fixations <- function(obj,rois='all',append=FALSE){

  allrois <- lapply(obj$rois,function(x) x$roi)
  allnames <- unlist(lapply(obj$rois, function(x) x$name))

  allnames <- paste0('roi_',allnames)

  hits <- data.frame(fixation_key = obj$fixations$fixation_key)

  names(hits)[-1] <- allnames

  for(i in 1:length(allrois)){
   hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$fixations$gavx,obj$fixations$gavy,allrois[[i]]))
  }

  hits <- hits[,c('fixation_key',allnames)]

  if(!append)
    tmp <- obj$fixations[c('fixation_key','ID','eyetrial','sttime','entime','gavx','gavy')]
  else
    tmp <- obj$fixations

  obj$fixations <- cbind(tmp,dplyr::select(hits,-fixation_key))

  return(obj)
}

calcHits_saccades <- function(obj,rois='all',append=FALSE){

  allrois <- lapply(obj$rois,function(x) x$roi)
  allnames <- unlist(lapply(obj$rois, function(x) x$name))

  allnames <- c(paste0('roi_start_',allnames), paste0('roi_end_',allnames))

  hits <- data.frame(saccade_key = obj$saccades$saccade_key)
  names(hits)[-1] <- allnames

  for(i in 1:length(allrois)){
    hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$saccades$gstx,obj$saccades$gsty,allrois[[i]]))
  }

  for(i in (length(allrois)+1):(length(allrois)*2)){
    hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$saccades$genx,obj$saccades$geny,allrois[[i-length(allrois)]]))
  }

  hits <- hits[,c('saccade_key',allnames)]

  if(!append)
    tmp <- obj$saccades[c('saccade_key','ID','eyetrial','sttime','entime','gstx','gsty','genx','geny')]
  else
    tmp <- obj$saccades

  obj$saccades <- cbind(tmp,dplyr::select(hits,-saccade_key))

  return(obj)
}

mapROIs <- function(obj,names,indicators=NULL,cleanup=FALSE){

  newnames <- names
  startnames <- newnames

  roinames <- unlist(lapply(obj$rois, function(x) x$name))

  for(eyedata in c('fixations','saccades')){

    if(eyedata=='fixations'){
      keyvar='fixation_key'
      roivars <- paste0('roi_',roinames)
      newnames <- paste0(startnames,'_hit')
    }
    else if(eyedata=='saccades'){
      keyvar='saccade_key'
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

epoch_fixations <- function(obj,roi,start=0,end=700,binwidth=25,event=NULL,type=NULL){

  startvar <- 'sttime'
  endvar <- 'entime'

#~   if(is.numeric(roi))
#~     hitvar <- paste0('roi_',roi)
#~   else
    hitvar <- paste0(roi,'_hit')

  prefix = 't'
  firstbin <- (ceiling(start/binwidth))
  numbins = ((end - start)/binwidth)+1
  lastbin <- firstbin + numbins-1

  if(is.null(event)){
    obj$header$default_event = obj$header$starttime
    event = 'default_event'
  }

  eventdata <- obj$header[,c('ID','eyetrial','starttime',event)]

  df <- obj$fixations[,c('ID','eyetrial','fixation_key',startvar,endvar,hitvar)]
  df <- merge(df,eventdata,by=c('ID','eyetrial'),all.x=T)
  df <- dplyr::arrange(df,fixation_key)

  #convert to "trial time"
  df[event] <- df[,event] - df$starttime
  df[startvar] <- df[,startvar] - df$starttime
  df[endvar] <- df[,endvar] - df$starttime

  #bin the data based on saccade start time and bin width (in ms)
  df$bin_start <- ceiling((df[,startvar] - df[,event])/binwidth)

  if(type=='cumulative'){
    df$bin_end <- NA
    df$bin_end[df$bin_start>=firstbin & df$bin_start<=lastbin] <- lastbin
  }
  else
    df$bin_end <- ceiling((df[,endvar] - df[,event])/binwidth)

  #throw out saccades before stim onset and after numbins
  df <- subset(df,bin_start>=firstbin & bin_end<=lastbin)

  #convert to "trial time"
#~   df[startvar] <- df[,startvar] - df$starttime
#~   df[endvar] <- df[,endvar] - df$starttime

  if(start<0){
    start_adj = df$bin_start + abs(firstbin)+1
    end_adj = df$bin_end + abs(firstbin)+1
    intervals <- apply(cbind(start_adj,end_adj),1,function(x) x[1]:x[2])
  }
  else {
    #create intervals (start_bin:end_bin)
    intervals <- apply(df[c('bin_start','bin_end')],1,function(x) x[1]:x[2])
  }

  #create a vector for each saccade (numbins wide). Fill with NA, then put 1's wherever saccade occurred (start_bin:end_bin)
  f2 <- function(x){
    vec <- rep(NA,numbins)
    vec[x] <-1
    vec
  }

  #this creates a list
  vectors <- lapply(intervals,f2)

  #stack into a matrix
  vectors <- do.call(rbind,vectors)

  #now, multiply each column by our vector indicating whether that saccade "hit" our item of interest
  #(will produce 0's where saccade occurred, but missed, 1 where it occured and missed, and NA where it didn't occur)
  vectors <- sweep(vectors,MARGIN=1,df[,hitvar],'*')

  #save the result in our data frame
  vars = gsub('-','_',paste0(prefix,firstbin:lastbin))
  df[vars] <- vectors

  #merge with our behavioral data
  df$roi <- roi
  df$epoch_start <- start
  df$epoch_end <- end
  df$binwidth <- binwidth
  df <- df[,c('ID','eyetrial','fixation_key','roi','epoch_start','epoch_end','binwidth',startvar,endvar,'bin_start','bin_end',event,hitvar,vars)]

  df <- merge(obj$beh[c('ID','eyetrial',obj$indexvars)],df,by=c('ID','eyetrial'),all.y=T)
  df <- dplyr::arrange(df,fixation_key)

  if(is.null(obj$epochs$fixations))
    obj$epochs$fixations <- list()

  if(is.null(obj$epochs$fixations[[event]]))
    obj$epochs$fixations[[event]] <- list()

  obj$epochs$fixations[[event]][[roi]] <- df

  # force garbage collection:
  gc()

  return(obj)
}

do_agg_fixations <- function(obj,event,roi,groupvars=c(),level='group',shape='long',condition="NULL"){

  prefix <- 't'

  # df <- obj$epochs$fixations[[event]][[roi]]

  df <- eyemerge(obj,'epoched_fixations',behdata=groupvars,event=event,roi=roi,condition=condition,condition.str=T)

  # remove unneccesary columns to save memory
  # df <- dplyr::select(df, c(-fixation_key, -sttime, -entime, -bin_start, -bin_end)) #, -epoch_start, -epoch_end, -binwidth))

  #remove specified filters: filter[1] == column name, filter[2] == condition (e.g., preposition == 'above')
  # if(!is.null(filter)){
    # var = filter[1]
    # val = filter[2]
    # filter_criteria <- lazyeval::interp(~ which_column == val, which_column = as.name(var))
    # df <- dplyr::filter_(df,filter_criteria)

    # df <-df[eval(parse(text=paste0("df$",filter[1]))) == filter[2], ]
  # if(condition!="NULL"){
  #   condition_call <- parse(text=condition)
  #   r <- eval(condition_call,df, parent.frame())
  #   df <- df[r, ]
  # }
  # # }



  epoch_start <- df$epoch_start[1]
  epoch_end <- df$epoch_end[1]
  binwidth <- df$binwidth[1]

  binnames <- names(df)[grepl(paste0('^',prefix,'[_0-9]'),names(df))]
  #reshape the data (turn our bins into rows)
  # na.rm saves memory
  df <- tidyr::gather_(df,'bin','val',binnames, na.rm = TRUE)

 #aggregate by trial(max)
  #this is super weird, just to make things play nice with dplyr
  varnames = sapply(c(obj$idvar,obj$indexvars,groupvars,'bin'), . %>% {as.formula(paste0('~', .))})
  df <- dplyr::group_by_(df,.dots=varnames) %>%
    dplyr::summarise(val = max(val,na.rm=T))

 if(level!='trial'){
    #aggregate by subject (mean)
    varnames2 = sapply(c(obj$idvar,groupvars,'bin'), . %>% {as.formula(paste0('~', .))})
    df <- dplyr::group_by_(dplyr::ungroup(df),.dots=varnames2) %>%
      dplyr::summarise(val=mean(val,na.rm=T))
  }
  if(level=='group'){
    #aggregate across subjects
    varnames3 = sapply(c(groupvars,'bin'), . %>% {as.formula(paste0('~', .))})
    df <-  dplyr::group_by_(dplyr::ungroup(df),.dots = varnames3) %>%
      dplyr::summarise(val=mean(val,na.rm=T))
  }

  df$roi <- roi
  df$epoch_start <- epoch_start
  df$epoch_end <- epoch_end
  df$binwidth <- binwidth

  if(shape=='wide'){
    df <- df %>%
      tidyr::spread(bin,val)
  }

  return(dplyr::ungroup(df))
}

aggregate_fixation_timeseries <- function(obj,event,rois,groupvars=c(),level='group',shape='long',type='probability',condition=NULL,condition.str=FALSE){

  agg <- data.frame()

  if(!condition.str)
    condition <- deparse(substitute(condition))

  if((type != 'probability') && length(rois==2)){

    if(any(groupvars == 'roi')){
      stop('Using "roi" as a grouping variable does not work if you want to plot a contrast score.')
    }

    aggID_one <- aggregate_fixation_timeseries(obj,event=event,roi=rois[1],groupvars = groupvars,shape='long',level='ID',condition=condition,condition.str=T)
    aggID_two <- aggregate_fixation_timeseries(obj,event=event,roi=rois[2],groupvars = groupvars,shape='long',level='ID',condition=condition, condition.str=T)

    # force garbage collection:
    gc()

    agg <- dplyr::inner_join(aggID_one,aggID_two,by=c('ID',groupvars,'bin'))

    if( ( (!all(agg$binwidth.x==agg$binwidth.y)) || (!all(agg$epoch_start.x==agg$epoch_start.y)) || (!all(agg$epoch_end.x==agg$epoch_end.y)) ))
      stop('Epochs/binning is different for the different ROIs. Data will not match up')

    agg <- dplyr::select(agg, -epoch_start.y, -epoch_end.y, -binwidth.y)

    agg <- dplyr::rename(agg,
                           epoch_start = epoch_start.x,
                           epoch_end = epoch_end.x,
                           binwidth = binwidth.x)

    if(type == 'difference'){
      agg$val <- agg$val.x - agg$val.y
    } else if (type == 'logRatio'){
      agg$val <- log((agg$val.x+1.0) / (agg$val.y+1.0))
    } else if(type == 'proportion'){
      agg$val <- agg$val.x / agg$val.y
    }

    #select only the variables we want
    agg <- dplyr::select_(agg, .dots= c('ID',groupvars,'bin','val','epoch_start','epoch_end','binwidth'))

    if(level=='group'){
      #group_by using string variable names
      dots <- lapply(c('bin','epoch_start','epoch_end','binwidth',groupvars), as.symbol)
      agg <- dplyr::group_by_(agg,.dots=dots) %>%
        dplyr::summarise(val = mean(val,na.rm=T))
    }

    agg$roi <- type

  }
  else{
    for(r in rois){
      aggtmp <- do_agg_fixations(obj,event=event,
                                              roi=r,
                                              groupvars = groupvars,
                                              shape=shape,
                                              level=level,
                                              condition=condition)


      if(nrow(agg)==0)
        agg <- aggtmp
      else
        agg <- rbind(agg,aggtmp)
    }
  }

  # force garbage collection:
  gc()

  return(agg)
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
