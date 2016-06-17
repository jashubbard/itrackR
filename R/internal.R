rois2df <- function(obj,which='all'){

  roidfs <- list()

  for(i in 1:length(obj$rois)){

    eltmp <- as.data.frame(spatstat::vertices(obj$rois[[i]]$roi))

    eltmp$name <- obj$rois[[i]]$name
    #as numeric because spatstat::centroid.owin (used for the center of polygons) seems to not return numerics
    eltmp$xcenter <- as.numeric(obj$rois[[i]]$center[1])
    eltmp$ycenter <- as.numeric(obj$rois[[i]]$center[2])
    eltmp$xradius <- obj$rois[[i]]$xradius
    eltmp$yradius <- obj$rois[[i]]$yradius
    eltmp$radius <- obj$rois[[i]]$radius
    eltmp$shape <- obj$rois[[i]]$shape

    roidfs[[i]]<- eltmp

  }

  df <- do.call('rbind',roidfs)
  return(df)


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


check_for_samples <- function(obj)
{

  #quickly check if all sample files exist. If so, just stop
  if(all(unlist(lapply(obj$samples,file.exists))) && length(obj$samples)==length(obj$edfs))
    return(obj)
  else{

    print('At least some sample data not loaded yet. Loading now (this make take a while)...')
    obj <- load_samples(obj)

    return(obj)

  }

}


timeshift <- function(df,baseline,cols){
  #adjusts the named columns by some baseline
  #baseline can be scalar or a vector
  #useful for computing relative time from start of experiment, or trial-by-trial time.

  df[cols] <- df[,cols] - baseline

  return(df)







}

