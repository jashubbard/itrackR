
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

  allrois <- lapply(obj$rois,function(x) x$roi)
  allnames <- unlist(lapply(obj$rois, function(x) x$name))


  if(is.numeric(allnames))
    allnames <- paste0('roi_',allnames)

  hits <- data.frame(fixation_key = obj$fixations$fixation_key)
  names(hits)[-1] <- allnames

  for(i in 1:length(allrois)){


   hits[allnames[[i]]] <- as.numeric(spatstat::inside.owin(obj$fixations$gavx,obj$fixations$gavy,allrois[[i]]))

  }

  hits <- hits[,c('fixation_key',paste0('roi_',1:length(obj$rois)))]

  if(!append)
    tmp <- obj$fixations[c('fixation_key','ID','eyetrial','sttime','entime','gavx','gavy')]
  else
    tmp <- obj$fixations

  obj$fixations <- cbind(tmp,dplyr::select(hits,-fixation_key))

  return(obj)

}


mapROIs <- function(obj,names,indicators=NULL){

  names <- paste0(names,'_hit')

  #if we're re-running this, first remove the existing columns
  if(any(names(obj$fixations) %in% names))
    obj$fixations <- obj$fixations[ , -which(names(obj$fixations) %in% names)]


  tmp1 <- obj$beh[,c('ID','eyetrial',indicators)] #behavioral data with trial-by-trial indicators
  tmp2 <- obj$fixations[,c('ID','eyetrial','fixation_key',paste0('roi_',1:length(obj$rois)))] #fixation data with roi hits
  tmp3 <- merge(tmp1,tmp2,by=c('ID','eyetrial'),all.y=T) #merge them together so we have our indicators for fixation-by-fixation
  tmp3 <- dplyr::arrange(tmp3,fixation_key) #put in the original order
  hits <- tmp3[,c(paste0('roi_',1:length(obj$rois)))] #grab only the hit data (since we will use the indicator to index from this)

  #funny hack to handle situations where indicator is NA.
  #Create a column at the end that is all NAs, so when indexed it fills in that value
  hits$missing <- NA
  hits <- as.matrix(hits)

  for(i in 1:length(indicators)){

    tmp3[names[i]] <- NA

    #wherever indicator is NA, change to #rois+1, so it indexes the "missing" column created above
    tmp3[,indicators[i]][is.na(tmp3[,indicators[i]])] <- ncol(hits)

    idx <- sub2ind(hits,1:nrow(tmp3),tmp3[,indicators[i]]) #get the single number indices based on row and column numbers (Indicator corresponds to column number)
    tmp3[names[i]] <- hits[idx] #grab the actual data based on the single number indices. Yields a single column of hits and misses

    tmp3[,indicators[i]][tmp3[,indicators[i]]==ncol(hits)] <- NA #fill the missing data back with NAs
  }

  #in case we're re-running this, don't create duplicate columns
  if(any(names(obj$fixations) %in% indicators))
      obj$fixations <- obj$fixations[ , -which(names(obj$fixations) %in% indicators)]


  #merge our existing fixation data with our new mapped data
  obj$fixations <- merge(obj$fixations,tmp3[,c('fixation_key',indicators,names)],by='fixation_key',all.x=T)

  return(obj)

}

sub2ind <- function(mat, r,c){

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

epoch_fixations <- function(obj,roi,start=0,end=700,binwidth=25,event='STIMONSET'){

  startvar <- 'sttime'
  endvar <- 'entime'
  hitvar <- paste0(roi,'_hit')
  prefix = 't'
  numbins = (end - start)/binwidth

  eventdata <- obj$header[,c('ID','eyetrial','starttime',event)]

  df <- obj$fixations[,c('ID','eyetrial','fixation_key',startvar,endvar,hitvar)]
  df <- merge(df,eventdata,by=c('ID','eyetrial'),all.x=T)
  df <- dplyr::arrange(df,fixation_key)
  df[startvar] <- df[,startvar] - df$starttime
  df[endvar] <- df[,endvar] - df$starttime
  df[event] <- df[,event] - df$starttime

  #bin the data based on saccade start time and bin width (in ms)
  df$bin_start <- ceiling((df[,startvar] - df[,event])/binwidth)
  df$bin_end <- ceiling((df[,endvar] - df[,event])/binwidth)
  #throw out saccades before stim onset and after numbins
  df <- subset(df,bin_start>start & bin_end<=numbins)


  #create intervals (start_bin:end_bin)
  intervals<- apply(df[c('bin_start','bin_end')],1,function(x) x[1]:x[2])

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
  vars = paste0(prefix,1:numbins)
  df[vars] <- vectors

  #merge with our behavioral data
  df$roi <- roi
  df$epoch_start <- start
  df$epoch_end <- end
  df$binwidth <- binwidth
  df <- df[,c('ID','eyetrial','fixation_key','roi','epoch_start','epoch_end','binwidth',startvar,endvar,hitvar,vars)]

  df <- merge(obj$beh,df,by=c('ID','eyetrial'),all.y=T)
  df <- dplyr::arrange(df,fixation_key)


  if(is.null(obj$epochs$fixations))
    obj$epochs$fixations <- list()

  if(is.null(obj$epochs$fixations[[event]]))
    obj$epochs$fixations[[event]] <- list()


  obj$epochs$fixations[[event]][[roi]] <- df

  return(obj)
}


do_agg_fixations <- function(obj,event,roi,groupvars,level='group',shape='long'){
  idvar <- 'ID'
  prefix <- 't'

  #reshape the data (turn our bins into rows)
  df <- obj$epochs$fixations[[event]][[roi]]

  epoch_start <- df$epoch_start[1]
  epoch_end <- df$epoch_end[1]
  binwidth <- df$binwidth[1]

  binnames <- names(df)[grepl(paste0(prefix,'[0-9]'),names(df))]
  df <- tidyr::gather_(df,'bin','val',binnames)

  #aggregate by trial(max)
  #this is super weird, just to make things play nice with dplyr
  varnames = sapply(c(idvar,obj$indexvars,groupvars,'bin'), . %>% {as.formula(paste0('~', .))})

  df <- dplyr::group_by_(df,.dots=varnames) %>%
    dplyr::summarise(val = max(val,na.rm=T))

  if(level!='trial' || level=='subject' || level=='ID'){
    #aggregate by subject (mean or median)
    varnames2 = sapply(c(idvar,groupvars,'bin'), . %>% {as.formula(paste0('~', .))})

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

  return(df)
}


aggregate_fixation_timeseries <- function(obj,event,rois,groupvars,level='group',shape='long',difference=FALSE){

  agg <- data.frame()

  if(difference && length(rois==2)){
    aggID_one <- aggregate_fixation_timeseries(obj,event='STIMONSET',roi=rois[1],groupvars = c('Task','Conflict'),shape='long',level='ID')
    aggID_two <- aggregate_fixation_timeseries(obj,event='STIMONSET',roi=rois[2],groupvars = c('Task','Conflict'),shape='long',level='ID')

    aggID <- merge(aggID_one,aggID_two,by=c('ID','Task','Conflict','bin'))

    if( ( (!all(aggID$binwidth.x==aggID$binwidth.y)) || (!all(aggID$epoch_start.x==aggID$epoch_start.y)) || (!all(aggID$epoch_end.x==aggID$epoch_end.y)) ))
      stop('Epochs/binning is different for the different ROIs. Data will not match up')

    aggID$val <- aggID$val.x - aggID$val.y

    aggID <- dplyr::rename(aggID,
                           epoch_start = epoch_start.x,
                           epoch_end = epoch_end.x,
                           binwidth = binwidth.x)

    aggID <- dplyr::select(aggID,ID,Task,Conflict,bin,val,epoch_start,epoch_end,binwidth)


    if(level=='group'){
      agg <- dplyr::group_by(aggID,bin,Task,Conflict,epoch_start,epoch_end,binwidth) %>%
        dplyr::summarise(val = mean(val,na.rm=T))
    }
    else
      agg <- aggID

    agg$roi <- 'difference'

  }
  else{

  for(r in rois){
    aggtmp <- do_agg_fixations(obj,event=event,
                                            roi=r,
                                            groupvars = groupvars,
                                            shape=shape,
                                            level=level)

    agg <- rbind(agg,aggtmp)
  }

  }

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

