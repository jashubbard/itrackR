# epoch_samples <- function(obj,event,field='pa',epoch=c(-100,100))
# {
#
#   obj <- check_for_samples(obj) #if we haven't loaded sample data yet, this will do it
#
#   # allepochs <- list()
#
#   if(!(event %in% names(obj$header)))
#     stop("Time-locking event not found. Have you run find_messages yet?")
#
#   for(i in 1:length(obj$edfs))
#   {
#
#     #     justfile <- sub("^([^.]*).*", "\\1", basename(obj$edfs[[i]]))
#     #     id <- as.numeric(gsub("([0-9]*).*","\\1",justfile))
#
#     id <- obj$subs[[i]]
#
#     header <- subset(obj$header,ID==id)
#     events <- header[[event]][!is.na(header[[event]])]
#     trials <- header$eyetrial[!is.na(header[[event]])]
#     samples <- read_saved_samples(obj$samples[[i]],ID=obj$subs[i],cols=c('ID','time','eyetrial',field))
#     # samples <- readRDS(obj$samples[[i]])
#
#     thisepoch <- edfR::epoch.samples(events,as.data.frame(samples),sample.field=field,epoch=epoch,eyetrial=T)
#
#     thisepoch$ID <- rep(id,length(trials))
#     thisepoch$event <- event
#     thisepoch$eyetrial <- trials
#     obj$epochs$samples[[event]][[i]] <- thisepoch
#
#     rm(samples)
#     gc()
#
#   }
#
#   # obj$epochs[[event]] <- allepochs
#
#   #record what we did
#   obj$history$step <- obj$history$step +1
#
#   funcall <- list()
#   funcall$name <- 'epoch_samples'
#   funcall$step <- obj$history$step
#   funcall$event <- event
#   funcall$field <- field
#   funcall$epoch <- epoch
#
#   if(is.null(obj$history$epoch_samples))
#     obj$history_epoch_samples <- funcall
#   else
#     obj$history$epoch_samples <- c(obj$history$epoch_samples,funcall)
#
#   return(obj)
#
# }

load_samples <- function(obj,outdir=NULL, force=F,parallel=TRUE, ncores = 2){

  allsamps <- list()

  #if we're not forcing and all the files are there, then just go back
  if(!force && all(unlist(lapply(obj$samples,file.exists))))
    return(obj)

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
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    #run with dopar
    allsamps <- foreach::foreach(edf=iterators::iter(1:length(obj$edfs)),.packages=c('data.table')) %dopar%
    {
      fname <- load_sample_file(obj,edf,force=force)
      fname
    }
    #stop cluster when done
    parallel::stopCluster(cl)
  }
  else{

    #serial method - use lapply
    allsamps <- lapply(1:length(obj$edfs), function(x) load_sample_file(obj,x,force=force))

  }


  obj$samples <- allsamps
  print("complete")
  return(obj)

}


#write each sample file as a separate database
load_sample_file <- function(obj,i,force=F){

  #create the directory if it doesn't exist
  dir.create(path.expand(obj$sample.dir), showWarnings = FALSE)
  #full path to file (need for loading .edf)
  dbname <- file.path(path.expand(obj$sample.dir),paste0(obj$subs[i],'_samples.sqlite'))
  #relative path to file (more convenient)
  dbname_relative <- file.path(obj$sample.dir,paste0(obj$subs[i],'_samples.sqlite'))

  #don't reload if it's already done and we don't set force=T
  if(!force && file.exists(dbname))
    return(dbname_relative)

  #load the samples
  samps <- edfR::edf.samples(obj$edfs[i],trials=T,eventmask=T)

  #if it's not binocular, don't keep all the L/R data
  if(!obj$binocular){
    if(all(is.na(samps$paL))){
      samps[,(c("paL","gxL","gyL")) := NULL] #in-place delete of column
      data.table::setnames(samps,c("paR","gxR","gyR"),c("pa","gx","gy"))

    } else{

      samps[,(c("paR","gxR","gyR")) :=NULL]
      data.table::setnames(samps,c("paL","gxL","gyL"),c("pa","gx","gy"))

    }
  }

  #time-shift so samples start at 1
  baseline <- samps[1,time]-1
  samps[,time := time-baseline]

  samps$ID <- obj$subs[i]

  #create a separate database file for each subject in sample.dir
  #one big db would be nice, but can't do in parallel
  db <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=dbname)
  RSQLite::dbGetQuery(db,'PRAGMA page_size=4096') #increase the page size for increased performance (but a bit more memory)

  #write all the data into a table called samples
  RSQLite::dbWriteTable(db,'samples',samps,overwrite=T,append=F,row.names=F)

  #create an index for time and for trial. Makes searching and joining faster
  RSQLite:: dbGetQuery(db,"CREATE INDEX subtime ON samples (ID,time)")
  RSQLite::dbGetQuery(db,'CREATE INDEX person ON samples (ID,eyetrial)')
  RSQLite::dbDisconnect(db)


  rm(samps) #clear memory
  gc()

  print(sprintf('Finished with samples from file %s',obj$edfs[i]))

  return(dbname_relative)
}


# load_sample_file_old <- function(obj,edf,force=FALSE,db=FALSE){
#
#    id <- edf2id(edf)
#
#    if(!db)
#     fname <- file.path(obj$sample.dir,paste0(id,'_samp.rds'))
#    else
#      fname <- file.path(obj$sample.dir,'samples.sqlite')
#
#   if(!db && file.exists(fname) && !force)
#     samps <- readRDS(fname)
#   else{
#     samps <- edfR::edf.samples(edf,trials=T,eventmask=T)
#
#
#     #if it's not binocular, don't keep all the L/R data
#     if(!obj$binocular){
#       if(all(is.na(samps$paL))){
#         samps[,(c("paL","gxL","gyL")) := NULL] #in-place delete of column
#         data.table::setnames(samps,c("paR","gxR","gyR"),c("pa","gx","gy"))
#
#       } else{
#
#         samps[,(c("paR","gxR","gyR")) :=NULL]
#         data.table::setnames(samps,c("paL","gxL","gyL"),c("pa","gx","gy"))
#
#       }
#     }
#
#     #time-shift
#     baseline <- samps[1,time]-1
#     samps[,time := time-baseline]
#
#     if(!db)
#        saveRDS(samps,fname,compress = T)
#
#     else{
#
#       db.pupil = RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname)
#
#       overwrite = F
#
#       if(file.exists(fname) && !force){
#         append=T
#       }
#       else{
#         append = F
#         overwrite = T
#       }
#
#       samps$ID <- id
#
#       #write data to database. Append if file already exists
#       RSQLite::dbWriteTable(conn=db.pupil, name='SAMPLES',samps,overwrite=overwrite,append=append,row.names=F)
#       RSQLite::dbDisconnect(db.pupil)
#
#     }
#
#   }
#
#   return(fname)
# }



interpolate.gaps <- function(y,gaps)
{

  ## Locate gaps, blink starts and blink ends
  # A blink starts either when the blink variable goes from 0 to 1
  # or if its first value is 1. Similarly, a blink ends when the blink
  # variable goes from 1 to 0 or when the last value is 1.

  y.na <- is.na(y)
  if (all(!y.na)) return(y) # if no missing data, return y
  blink.start <- c(which.max(gaps), which(diff(gaps)==1) + 1)
  blink.start <- unique(blink.start) # remove eventual duplicates
  blink.end <- c(which(diff(gaps)==-1), max(which(gaps==1)))
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




# remove_blinks <- function(obj, interpolate=FALSE)
# {
#
#   obj <- check_for_samples(obj)
#
#   for(i in 1:length(obj$edfs))
#   {
#
#     samps <- read_saved_samples(obj$samples[[i]],ID=obj$subs[i],cols=c('ID','time','gx','gy','pa'))
#     # samps<- readRDS(obj$samples[[i]])
#
#     #code bad samples as "blinks"
#     samps$blink[samps$gx==1e08 | samps$gy==1e08] <- 1
#
#     blinks <- as.logical(samps$blink)
#     samps$pa[blinks] <-NA
#     samps$gx[blinks] <-NA
#     samps$gy[blinks] <-NA
#
#     if(interpolate){
#
#       samps[, pa := interpolate.gaps(samps$pa,samps$blink)]
#       samps[, gx := interpolate.gaps(samps$gx,samps$blink)]
#       samps[, gy := interpolate.gaps(samps$gy,samps$blink)]
#     }
#
#
#
#     saveRDS(samps,obj$samples[[i]],compress = T)
#     rm(samps)
#
#   }
#
#   #record what we did
#   obj$history$step <- obj$history$step + 1
#
#   funcall <- list()
#   funcall$name <- 'remove_blinks'
#   funcall$step <- obj$history$step
#   funcall$interpolate <- interpolate
#
#
#   return(obj)
# }



# get_all_epochs <- function(obj,epochname,baseline=NULL,baseline.method='percent',shape='wide',beh='all')
# {
#
#   window_start <- obj$epochs$samples[[epochname]][[1]]$epoch.window[1]
#   window_end <- obj$epochs$samples[[epochname]][[1]]$epoch.window[2]
#
#   epochs <- obj$epochs$samples[[epochname]]
#   epochs <- lapply(epochs,function(x) cbind(x$ID,x$eyetrial,x$epochs))
#   epochs <- do.call(rbind,epochs)
#
#   if(!is.null(baseline))
#   {
#     epochs_b <- baseline_epochs(epochs[,-c(1,2)],baseline=baseline,method=baseline.method)
#     epochs <- cbind(epochs[,1:2],epochs_b)
#   }
#
#   timenames <- paste0('t',seq(window_start,window_end-1))
#   timenames <- gsub('-','_',timenames)
#   colnames(epochs) <- c('ID','eyetrial',timenames)
#
#   epochs <- as.data.frame(epochs)
#
#   if(shape=='long'){
#
#     epochs <- tidyr::gather_(epochs,'timepoint','value',timenames)
#     epochs$timepoint <- gsub('t','',epochs$timepoint)
#     epochs$timepoint <- as.numeric(gsub('_','-',epochs$timepoint))
#   }
#
#   if(nrow(obj$beh)==0)
#     warning('No behavioral data was found. If you want to include it, run add_behdata')
#
#   if(!is.null(beh) && nrow(obj$beh)>0)
#   {
#
#     if(beh[1]=='all' || (is.logical(beh) && beh==TRUE))
#       behnames <- names(obj$beh)
#     else
#       behnames <- unique(c('ID','eyetrial',beh))
#
#     # behonly <- dplyr::select_(obj$beh,.dots=behnames)
#     behonly <- obj$beh[,behnames]
#     epochs <- dplyr::right_join(behonly,epochs,by=c('ID','eyetrial'))
#
#     if(shape=='long')
#       epochs <- dplyr::arrange(epochs,ID,eyetrial,timepoint)
#     else
#       epochs <- dplyr::arrange(epochs,ID,eyetrial)
#
#   }
#
#   return(epochs)
# }


# baseline_epochs <- function(epochs,baseline=c(1,100),method='percent'){
#
#
#   mn <- rowMeans(epochs[,baseline[1]:baseline[2]],na.rm=T)
#
#   epoch_b <- sweep(epochs,1,mn,FUN='-')
#   epoch_b <- sweep(epoch_b,1,mn,FUN='/')
#
#   return(epoch_b)
# }



# read_saved_samples <- function(fname,ID=NULL,cols = NULL){
#
#   if(tolower(tools::file_ext(fname))=='rds')
#     data <- readRDS(fname)
#
#   if(tolower(tools::file_ext(fname))=='sqlite'){
#
#     my_db <- dplyr::src_sqlite(fname, create = F)
#
#     data <- tbl(my_db, from = 'SAMPLES') %>%
#       filter_(.dots=c(paste0('ID==',ID)))
#
#     if(!is.null(cols))
#       data <- select_(data,.dots=cols)
#
#     data <- data.table::as.data.table(data)
#
#       rm(my_db)
#       gc()
#   }

#
#
#   return(data)
# }



get_sample_epochs <-function(obj,parallel=T,ncores=2,...){

if(length(obj$beh)==0)
  obj <- add_behdata(obj)

#create a database with behaivoral data
db2 <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=file.path(path.expand(obj$sample.dir),'beh.sqlite'))
RSQLite::dbWriteTable(db2,'beh',as.data.frame(obj$beh),row.names=F,overwrite=T)
RSQLite::dbGetQuery(db2,'CREATE INDEX person ON beh (ID,eyetrial)')
RSQLite::dbDisconnect(db2)

if(parallel){

  print('Getting sample epochs (in parallel)...')


  #set up cluster with maximum number of cores

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  #run with %dopar%
  epochs <- foreach::foreach(i=iterators::iter(1:length(obj$edfs)),.packages=c('data.table','itrackR'),.errorhandling = 'remove',.combine=dplyr::bind_rows) %dopar%
  {
    eps <- get_avg_epochs(obj,i,...)
    eps
  }

  parallel::stopCluster(cl)
}
else{

  #serial version - exactly the same, but with %do%
  epochs <- foreach::foreach(i=iterators::iter(1:length(obj$edfs)),.packages=c('data.table','itrackR'),.errorhandling = 'remove',.combine=dplyr::bind_rows) %do%
  {
    eps <- get_avg_epochs(obj,i,...)
    eps
  }
}

print('DONE!')


return(epochs)

}

get_avg_epochs <- function(obj,s,event='starttime', epoch = c(-100,100),factors=c('ID'), aggregate=T, baseline=NULL,cleanup=F){

  #get the subject ID
  subID <- obj$subs[s]

  #take the header of itrackr object, get only this subject's data
  #compute 2 variables, for the start and end time relative to the time-locking event
  tmp <- obj$header
  tmp <- dplyr::filter(tmp,ID==subID)
  tmp <- dplyr::filter_(tmp,paste0('!is.na(',event,')'))
  tmp <- dplyr::arrange(tmp,ID,eyetrial)
  tmp <- dplyr::mutate_(tmp,'tstart' = paste0(event,'+',epoch[1]),
                        'tend' = paste0(event,'+',epoch[2]))

  #for our queries, we need our factors as "Task, Conflict", not c("Task", "Conflict")
  facnames <- paste(factors,collapse = ', ')
  facnames_beh <- paste('beh.',factors,sep='',collapse=', ')  #to make "beh.Task, beh.Conflict"

  #connect to this subject's database
  dbname <- obj$samples[[s]]
  db <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=dbname)


  #loop through each row of the header, get the start and end time, extract the data we want
  for(i in 1:nrow(tmp)){

    if(i==1){
      #first row, create a temporary table called EPOCHS
      if(RSQLite::dbExistsTable(db,'EPOCHS'))
        RSQLite::dbGetQuery(db,'DROP TABLE EPOCHS')

      RSQLite::dbGetQuery(db,sprintf('CREATE TEMPORARY TABLE EPOCHS AS SELECT ID,time,eyetrial,pa, time - %d+1 AS timepoint, %d as epochnum FROM samples WHERE time BETWEEN %d AND %d AND ID=%d',tmp$tstart[i],i,tmp$tstart[i],tmp$tend[i],subID))
      RSQLite::dbBegin(db) #faster INSERT if we use BEGIN...COMMIT statement

    }else #otherwise, append to the already made EPOCHS table
      RSQLite::dbGetQuery(db,sprintf('INSERT INTO EPOCHS (ID,time,eyetrial,pa,timepoint,epochnum) SELECT ID,time,eyetrial,pa,time - %d+1 AS timepoint, %d as epochnum FROM samples WHERE time BETWEEN %d AND %d AND ID=%d',tmp$tstart[i],i, tmp$tstart[i],tmp$tend[i],subID))
  }

  #COMMIT when done inserting
  RSQLite::dbCommit(db)

  #if we want our data baselined
  if(!is.null(baseline)){

    #convert from relative-to-event to relative-to-epoch
    ep <- epoch[1]:epoch[2]
    bl <- c(which.max(baseline[1]==ep), which.max(baseline[2]==ep))

    if(!any(bl))
      stop('Problem with baseline period. This should be relative to your time-locking event, like your epoch')

    #create a new table, BASELINE that holds the baselined pupil data
    if(RSQLite::dbExistsTable(db,'BASELINE'))
      dbGetQuery(db,'DROP TABLE BASELINE')

    RSQLite::dbGetQuery(db,sprintf('CREATE TEMPORARY TABLE BASELINE AS SELECT *, (pa-baseline)/baseline as pupil_baseline FROM EPOCHS LEFT JOIN (SELECT ID,eyetrial,AVG(pa) as baseline FROM EPOCHS WHERE timepoint BETWEEN %d and %d GROUP BY epochnum) base ON base.ID=EPOCHS.ID AND base.eyetrial=EPOCHS.eyetrial',bl[1],bl[2]))
  }

  #for aggregating
  if(aggregate){
    #this holds everyone's behavioral data
    RSQLite::dbGetQuery(db,sprintf('ATTACH "%s" as beh', file.path(path.expand(obj$sample.dir),'beh.sqlite')))

    if(!is.null(baseline)){
      #if baselined, merge BASELINE with the behavioral data and compute the average by factors and timepoint
      res <- RSQLite::dbGetQuery(db,sprintf('SELECT ID,%s,timepoint, AVG(pupil_baseline) as pupil, AVG(pa) as pupil_raw FROM (SELECT BASELINE.ID,BASELINE.timepoint, BASELINE.eyetrial,BASELINE.pupil_baseline,BASELINE.pa,%s FROM BASELINE INNER JOIN beh ON beh.ID=BASELINE.ID AND beh.eyetrial=BASELINE.eyetrial) GROUP BY %s,timepoint ORDER BY %s',facnames,facnames_beh,facnames,facnames))

    }
    else{
      #otherwise, merge EPOCHS with the behavioral data and compute the average by factors and timepoint
      res <- RSQLite::dbGetQuery(db,sprintf('SELECT ID,%s,timepoint, AVG(pa) as pupil FROM (SELECT EPOCHS.ID,EPOCHS.timepoint, EPOCHS.eyetrial,EPOCHS.pa,%s FROM EPOCHS INNER JOIN beh ON beh.ID=EPOCHS.ID AND beh.eyetrial=EPOCHS.eyetrial) GROUP BY %s,timepoint ORDER BY %s',facnames,facnames_beh,facnames,facnames))
    }

    #don't need this anymore
    RSQLite::dbGetQuery(db,'DETACH beh')
  }
  else{
    #if not aggregating, give me all the data
    if(!is.null(baseline))
      res <- RSQLite::dbGetQuery(db,'SELECT * FROM BASELINE')
    else
      res <- RSQLite::dbGetQuery(db,'SELECT * FROM EPOCHS')
  }

  #remove the EPOCHS table
  if(cleanup){
    RSQLite::dbGetQuery(db,'DROP TABLE EPOCHS')
  }

  #disconnect from database
  RSQLite::dbDisconnect(db)

  #add the actuall time
  ep <- epoch[1]:epoch[2]
  res$epoch_time <- as.numeric(ep[res$timepoint])

  res$pupil <- as.numeric(res$pupil)

  if(!is.null(baseline)){
    res$pupil_raw <- as.numeric(res$pupil_raw)

    if(!aggregate)
      res$baseline <- as.numeric(res$baseline)

  }

  return(res)

}


remove_blinks <- function(obj,interpolate=F,parallel=T,ncores=2){

  if(parallel){

    print('Removing blinks (in parallel)...')


    #set up cluster with maximum number of cores

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)


    #run with dopar
    done <- foreach::foreach(dbname=iterators::iter(obj$samples),.packages = c('itrackR'), .errorhandling = 'remove') %dopar%
  {
    done <- blink_remove(dbname,interpolate=interpolate)
    done
  }

  parallel::stopCluster(cl)
}
else{

  print('Removing blinks...')

  #serial method - use lapply
  done <- lapply(obj$samples, function(x) blink_remove(x,interpolate=interpolate))
}

  return(obj)
}



blink_remove <- function(dbname,interpolate=FALSE){

  #connect to database file
  db <- DBI::dbConnect(RSQLite::SQLite(),dbname)

  #clear out blink or artifact periods in original data
  DBI::dbGetQuery(db,'UPDATE samples SET pa=NULL WHERE blink=1')

  if(interpolate){
    #pull out pupil data
    tmp <- DBI::dbGetQuery(db,'SELECT ID,time,pa,blink FROM samples')

    #remove data where blinks occur, then interpolate
    tmp$pa[tmp$blink==1] <- NA
    tmp$interp <- interpolate.gaps(tmp$pa,tmp$blink)

    #select out only the blink periods, and the ID, time, and interpolated data
    towrite <- dplyr::filter(tmp,blink==1)
    towrite <- dplyr::select(towrite,ID,time,interp)

    #add it to the database
    DBI::dbWriteTable(db,'pupil_interp',towrite,row.names=F,overwrite=T)
    #create index for speed
    DBI::dbGetQuery(db,'CREATE INDEX stime ON pupil_interp (ID,time)')

    #add our interpolated data to the blink periods
    DBI::dbGetQuery(db,'UPDATE samples SET pa = (SELECT interp FROM pupil_interp WHERE samples.ID=pupil_interp.ID AND samples.time=pupil_interp.time) WHERE blink=1')

    #thorw away the table
    DBI::dbRemoveTable(db,'pupil_interp')
  }

  return(DBI::dbDisconnect(db))

}



get_samples <- function(obj,IDs,fields=NULL,condition=NULL,parallel=F,ncores=2){


  idx <- which(obj$subs %in% IDs)

  if(parallel && length(idx)>1){

    print('Grabbing samples...')
    #set up cluster with maximum number of cores

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)

  results <- list()

  #run with dopar
  results <- foreach::foreach(dbname=iterators::iter(obj$samples[idx])) %dopar%
  {
    res <- load_sample_db(dbname,fields=fields,condition=condition)
    res
  }

  parallel::stopCluster(cl)
  results <- dplyr::bind_rows(results)
  }
  else
  {
    results <- load_sample_db(obj$samples[[idx]],fields=fields,condition=condition)

  }



  return(results)
}


load_sample_db <- function(dbname,fields,condition){

  db <- DBI::dbConnect(RSQLite::SQLite(),dbname=dbname)

  if(is.null(fields) && !is.null(condition))
    results <- DBI::dbReadTable(db,'samples')
  else{

    if(is.null(fields))
      sqlnames <- '*'
    else{
      fields <- fields[fields %in% RSQLite::dbListFields(db,'samples')]
      sqlnames <- paste(fields,collapse=', ')
    }

    query <- sprintf('SELECT %s FROM samples',sqlnames)

    if(!is.null(condition))
      query <- paste(query,sprintf('WHERE %s',condition))

    results <- DBI::dbGetQuery(db,query)
  }
  DBI::dbDisconnect(db)

  return(results)
}


check_interp <- function(obj,ID,start = NULL, end=NULL,window=5000){

  if(!is.null(start)){

    if(is.null(end))
      end <-start+window

    condition = sprintf('time>=%d AND time<=%d',start,end)

  }
  else
    condition=NULL

  samp <- get_samples(obj,IDs=ID,fields=c('pa','blink'),condition=condition)

  y <- samp$pa
  gaps <- samp$blink

  if(is.null(start) && is.null(condition)){
    warning('Time-locking to the first blink that I find...')
    start <- which.max(gaps==1)
  }
  else
    start <- 1

  if(is.null(end) && is.null(condition))
    end <- start+window
  else
    end <- nrow(samp)


  plot(y[start:end],type='l',col='black',lwd=1,xlab='timepoint',ylab='pupil diameter')
  par(new=T)

  y2 = y
  y2[gaps==0] <- NA

  lines(y2[start:end],col='red',lwd=3)


}
