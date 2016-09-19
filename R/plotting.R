# plot <- function(obj,zoom=FALSE,crosshairs=TRUE,rois=TRUE,which='all',names=FALSE){
#   UseMethod('plot')
#
# }

plot.itrackR <- function(obj,zoom=TRUE,crosshairs=TRUE,rois=TRUE,whichROIs='all',names=FALSE,IDs=c(),summarize=0,quick=F,condition=NULL){


  condition <- deparse(substitute(condition))

  if(condition !='NULL'){
    fixdata <- eyemerge(obj,condition=condition,condition.str = T)
  }
  else
    fixdata <- obj$fixations


  if (length(IDs) != 0)
    fixdata = dplyr::filter(fixdata,ID %in% IDs)

  if(quick)
    fixdata <- dplyr::distinct(dplyr::select_(fixdata,.dots=c(obj$idvar,'gavx','gavy')))

  if(rois && !is.null(obj$rois)){
    df <- rois2df(obj)

    if(whichROIs[1] !='all')
      df <- subset(df,name %in% whichROIs)

    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)

    if(summarize > 0){
      p <- p + ggplot2::stat_bin2d(data=fixdata,ggplot2::aes(x=gavx,y=gavy, size = ..count..), color='yellow',bins = summarize, geom = "point")
    } else {
      p <- p + ggplot2::geom_point(data=fixdata,ggplot2::aes(x=gavx,y=gavy),color='yellow',size=0.7)
    }

    if(names)
      p <- p + ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),aes(x=xcenter,y=ycenter,label=name),color='white')
  }
  else{
     if(summarize > 0){
      p <- ggplot2::ggplot(fixdata,ggplot2::aes(x=gavx,y=gavy)) + ggplot2::stat_bin2d(bins = summarize, ggplot2::aes(size = ..count..), color = 'yellow', geom = "point")
    } else {
      p <- ggplot2::ggplot(fixdata,ggplot2::aes(x=gavx,y=gavy)) + ggplot2::geom_point(color='yellow',size=0.7)
    }
  }

  if(zoom)
    p <- p + ggplot2::coord_cartesian(xlim=c(0,obj$resolution[1]), ylim=c(obj$resolution[2],0))
  else
    p <- p + ggplot2::coord_cartesian(xlim=c(0,max(fixdata$gavx)), ylim=c(max(fixdata$gavy),0))

  p <- p + ggplot2::scale_y_reverse() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank())

  if(crosshairs){
    p <- p +
      ggplot2::geom_hline(yintercept=round(obj$resolution[2]/2),color='red',size=0.3,linetype='dashed') +
      ggplot2::geom_vline(xintercept=round(obj$resolution[1]/2),color='red',size=0.3,linetype='dashed')
  }

  if(!summarize){
    p <- p + ggplot2::facet_wrap(~ID)
  }

  p
  return(p)
}

plot.timeseries <- function(obj,event=NULL,rois,lines,rows=NULL,cols=NULL,level='group',type='probability',condition=NULL){

  if(is.null(event)){
    event = 'default_event'
  }

  condition <- deparse(substitute(condition))

  agg <- aggregate_fixation_timeseries(obj,event=event,
                                       rois=rois,
                                       groupvars = c(lines,rows,cols),
                                       shape='long',
                                       level=level,
                                       type=type,
                                       condition=condition,
                                       condition.str = T)

  mainvars <- c('bin','val','roi','epoch_start','epoch_end','binwidth')
  othervars <- names(agg)[-which(names(agg) %in% mainvars)]

  tmp <- agg
  tmp$bin <- gsub('t','',tmp$bin)
  tmp$bin <- as.numeric(gsub('_','-',tmp$bin))
  tmp$bin <- (tmp$bin)*tmp$binwidth

  if(('roi' %in% lines) && (length(lines)>1))
    plotlines <- lines[-which('roi' %in% lines)]
  else
    plotlines <- lines


  tmp$lines <- interaction(agg[,plotlines])

  if(!is.null(rows)){
    tmp$rows <- interaction(agg[,rows])
    rowexp <- do.call('paste',c(as.list(rows),sep=' + '))
  }

  if(!is.null(cols)){
    tmp$cols <- interaction(agg[,cols])
    colexp <- do.call('paste',c(as.list(cols),sep=' + '))
  }



  plt <- ggplot2::ggplot(data=tmp,ggplot2::aes(x=bin,y=val,group=interaction(roi,lines),color=lines,shape=roi,linetype=roi))+
    ggplot2::geom_line() + #size=.75) +
    ggplot2::geom_point(size=1.8) +
    ggplot2::xlab('time (ms)') +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), alpha=0.5, linetype='dashed') +
    ggplot2::scale_color_discrete(name = toString(lines)) +
    ggplot2::scale_x_continuous() #breaks = seq(min(tmp$bin),max(tmp$bin),100))
    # using different than default breaks results in too much breaks (at least for me)

  if(type == 'difference'){
    plt <- plt + ggplot2::geom_hline(yintercept=0,linetype='dashed',alpha=0.7) +
            ggplot2::ylab('difference of fixations')
  } else if (type == 'logRatio'){
    plt <- plt + ggplot2::geom_hline(yintercept=0) +
            ggplot2::ylab('log gaze of fixations')
  } else if(type == 'proportion') {
    plot <- plt + ggplot2::ylab('proportion of fixations')
  } else if(type == 'probability') {
    plt <- plt + ggplot2::ylab('probability of fixations (%)')
  }

  if (!is.null(rows) | !is.null(cols)){

    if(is.null(rows))
      rowexp='.'
    if(is.null(cols))
      colexp='.'

    facetexp <- paste(rowexp,colexp,sep=" ~ ")

    plt <- plt + ggplot2::facet_grid(facetexp,labeller = ggplot2::label_both)
  }

  plt <- plt + ggplot2::theme(panel.background = ggplot2::element_blank(),
                              panel.border = ggplot2::element_blank(),
                              axis.line.x = ggplot2::element_line(color="black", size = .5),
                              axis.line.y = ggplot2::element_line(color="black", size = .5),
                              axis.text.x = ggplot2::element_text(size=14,vjust=-0.9),
                              axis.title.x = ggplot2::element_text(size=16,vjust=-1.5),
                              axis.text.y = ggplot2::element_text(size=14),
                              axis.title.y = ggplot2::element_text(size=16,vjust= -1.5),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())

    output = list()
    output$plt <- plt
    output$data <- tmp

    plt

  return(output)
}

plot_random_epochs <- function(epochs,n=100)
{

  if('timepoint' %in% names(epochs))
  {
    epochs <- dplyr::select(epochs,ID,eyetrial,timepoint,value)
    epochs <- tidyr::spread(epochs,timepoint,value)
  }


  epochs <- dplyr::select(epochs,matches('^t[_0-9]'))

  firstpoint <- gsub(colnames(epochs)[1],'t','')
  firstpoint <- as.numeric(gsub(firstpoint,'_','-'))

  lastpoint <- gsub(colnames(epochs)[ncol(epochs)],'t','')
  lastpoint <- as.numeric(gsub(lastpoint,'_','-'))

  matplot(t(as.matrix(epochs)),type='l',lty=1,axes=FALSE)

  axis(2)
  axis(1,at=seq(firstpoint,lastpoint,100),labels=seq(firstpoint,lastpoint,100)-(firstpoint))
  abline(v=firstpoint,col=4,lty=2,lwd=2)

}

plot.rois <- function(obj,which='all',crosshairs=T){

  df <- rois2df(obj)

  if(which[1] !='all')
    df <- subset(df,name %in% which)

  p <- ggplot2::ggplot() + ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)+
    ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),ggplot2::aes(x=xcenter,y=ycenter,label=name),color='white') +
    ggplot2::coord_cartesian(xlim=c(0,obj$resolution[1]), ylim = c(obj$resolution[2],0)) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  if(crosshairs){
    p <- p +
      ggplot2::geom_hline(yintercept=obj$resolution[2]/2,color='red',size=0.3,linetype='dashed') +
      ggplot2::geom_vline(xintercept=obj$resolution[1]/2,color='red',size=0.3,linetype='dashed')
  }

  p
  return(p)
}

plot.samples <- function(obj,ID,events=T,timestamp=NULL,showmean=T,bin=F,time.start=NULL,time.end=NULL){

  #load sample data if you haven't already
  obj <- check_for_samples(obj)

  #find the appropriate file that matches the ID you want
  sampfile <- unlist(obj$samples[grepl(ID,obj$samples)])

  #load in sample data
  samps <- readRDS(sampfile)
  samps <- as.data.frame(samps)
  samps$index <- 1:nrow(samps)

  timestart <- samps$time[1]

  if(!is.null(time.start))
    range_start <- timestart+time.start
  else
    range_start <- timestart


  if(!is.null(time.end))
    range_end <- timestart+time.end
  else
    range_end <- range_end <- max(samps$time)


  samps <- subset(samps,time>=range_start & time<=range_end)

  #for now, bin our sample data into 10 equal bins
  samps$bin <- findInterval(1:nrow(samps),seq(1,nrow(samps),nrow(samps)/10))


  plt <- ggplot2::ggplot()

  #downsample to every 100 samples
  downsamps <- samps[downsample(1:nrow(samps),100),]

  #limits for the y axes (5th precentile - 300 to maximum value +100)
  ylow <- quantile(samps$pa,.05, na.rm=T)-300
  yhigh <- max(samps$pa,na.rm=T)+100
  barheight <- (yhigh - ylow)/8


  if(showmean){
    #get the mean of all pupil data
    meanpupil <- mean(samps$pa,na.rm=T)
    #plot a horizontal dotted line
    plt <- plt + ggplot2::geom_hline(yintercept=meanpupil,alpha=0.5, linetype='longdash',size=0.5)
  }

  #draw little rectangles when certain events occur
  if(events){

    fixations <- obj$fixations[obj$fixations$ID==ID,]
    fixations$time <- fixations$sttime

    fixations <- subset(fixations,time>=range_start & time<=range_end)

    fixations <- dplyr::left_join(fixations,samps[c('time','index','bin')],by='time')
    fixations$ymin <- ylow
    fixations$ymax <- ylow+barheight

    fixations$xmin <- fixations$index
    fixations$xmax <- fixations$index + (fixations$entime - fixations$sttime)

    blinks <- obj$blinks[obj$blinks$ID==ID,]
    blinks$time <- blinks$sttime

    blinks <- subset(blinks,time>=range_start & time<=range_end)

    blinks <- dplyr::left_join(blinks,samps[c('time','index','bin')],by='time')
    blinks$ymin <- ylow+barheight+10
    blinks$ymax <- blinks$ymin+barheight
    blinks$xmin <- blinks$index
    blinks$xmax <- blinks$index + (blinks$entime - blinks$sttime)

    #do the actual plotting
    plt <- plt +
      ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,ymin=ymin,ymax=ymax),fill='red',alpha=0.6,data=blinks) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,ymin=ymin,ymax=ymax),fill='green',alpha=0.6,data=fixations)

  }

  #draw dotted lines that correspond to some timestamp
  if(!is.null(timestamp)){
    msgtime <- obj$messages[obj$messages$ID==ID & grepl(pattern = timestamp ,obj$messages$message),]

    if(nrow(msgtime)>0){
    msgtime$flag <- 1
    msgtime$time <- msgtime$sttime - timestart

    msgtime <- subset(msgtime,time>=range_start & time<=range_end)

    msgtime <- dplyr::left_join(msgtime,samps[c('time','index','bin')],by='time')
    msgtime <- subset(msgtime, !is.na(bin))
    plt <- plt + ggplot2::geom_vline(ggplot2::aes(xintercept=index),size=0.3,color='blue',linetype = "longdash",alpha=0.6,data=msgtime)
    }
    else
      warning('timestamp message not found')
  }

  xbreaks <- downsamps$index[seq(1,nrow(downsamps),100)]
  xlabels <- floor((downsamps$time[seq(1,nrow(downsamps),100)]-timestart)/1000)

  #draw the pupil data
    plt <- plt +
      ggplot2::geom_line(data=downsamps, ggplot2::aes(x=index,y=pa),size=1.0)  +
      ggplot2::theme_bw() +
      ggplot2::ylim(c(ylow,yhigh)) +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white'),
                              panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())

    if(bin)
      plt <- plt + ggplot2::facet_wrap(~bin,ncol=1,scales = 'free_x') + ggplot2::ylab('pupil size (arbitrary units)') + ggplot2::scale_x_continuous(name='time (seconds)',breaks=xbreaks,labels=xlabels)

    plt

    return(plt)
}

