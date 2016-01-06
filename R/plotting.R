plot <- function(obj,zoom=FALSE,crosshairs=TRUE,rois=TRUE,which='all',names=FALSE){
  UseMethod('plot')

}

plot.itrackR <- function(obj,zoom=FALSE,crosshairs=TRUE,rois=TRUE,which='all',names=FALSE){


  if(rois && !is.null(obj$rois)){
    df <- rois2df(obj)

    if(which !='all')
      df <- subset(df,name %in% which)


    p <- ggplot2::ggplot() + ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)+
      ggplot2::geom_point(data=obj$fixations,ggplot2::aes(x=gavx,y=gavy),color='yellow',size=0.7)

    if(names)
      p <- p + ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),aes(x=xcenter,y=ycenter,label=name),color='white')
  }
  else
    p <- ggplot2::ggplot(obj$fixations,ggplot2::aes(x=gavx,y=gavy)) + ggplot2::geom_point(color='yellow',size=0.7)


  if(zoom)
    p <- p+ggplot2::xlim(xlim=c(0,obj$resolution[1])) + ggplot2::ylim(ylim=c(obj$resolution[2],0))
  else
    p <- p + ggplot2::xlim(xlim=c(0,max(obj$fixations$gavx))) + ggplot2::ylim(ylim=c(max(obj$fixations$gavy),0))


  p <- p + ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank())


  if(crosshairs){
    p <- p +
      ggplot2::geom_hline(yintercept=round(z$resolution[2]/2),color='red',size=0.3,linetype='dashed') +
      ggplot2::geom_vline(xintercept=round(z$resolution[1]/2),color='red',size=0.3,linetype='dashed')
  }



  p <- p + ggplot2::facet_wrap(~ID)

  p
  return(p)
}


plot.timeseries <- function(obj,event,rois,lines,rows=NULL,cols=NULL,level='group',difference=FALSE){

  #   agg <- data.frame()
  #
  #   for(r in rois){
  #   aggtmp <- aggregate_fixation_timeseries(obj,event=event,
  #                                        roi=r,
  #                                        groupvars = c(lines,rows,cols),
  #                                        shape='long',
  #                                        level=level)
  #
  #   agg <- rbind(agg,aggtmp)
  #   }

  agg <- aggregate_fixation_timeseries(obj,event=event,
                                       rois=rois,
                                       groupvars = c(lines,rows,cols),
                                       shape='long',
                                       level=level,
                                       difference=difference)

  mainvars <- c('bin','val','roi','epoch_start','epoch_end','binwidth')
  othervars <- names(agg)[-which(names(agg) %in% mainvars)]

  tmp <- agg
  tmp$bin <- as.numeric(gsub('t','',tmp$bin))
  tmp$bin <- tmp$epoch_start + (tmp$bin-1)*tmp$binwidth
  tmp$lines <- interaction(agg[,lines])

  if(!is.null(rows)){
    tmp$rows <- interaction(agg[,rows])
    rowexp <- do.call('paste',c(as.list(rows),sep=' + '))
  }

  if(!is.null(cols)){
    tmp$cols <- interaction(agg[,cols])
    colexp <- do.call('paste',c(as.list(cols),sep=' + '))
  }


  plt <- ggplot2::ggplot(data=tmp,ggplot2::aes(x=bin,y=val,group=lines,color=lines))+
    ggplot2::geom_line(size=.75) +
    ggplot2::geom_point(size=3) + ggplot2::xlab('time (ms)') + ggplot2::ylab('probability of fixation (%)') +
    ggplot2::scale_color_discrete(name = toString(lines))


  if (!is.null(rows) | !is.null(cols)){

    if(is.null(rows))
      rowexp='.'
    if(is.null(cols))
      colexp='.'

    facetexp <- paste(rowexp,colexp,sep=" ~ ")


    plt <- plt + ggplot2::facet_grid(facetexp)
  }

  plt <- plt +  ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white'),
                               panel.border = ggplot2::element_blank(),
                               axis.line = ggplot2::element_line(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank())

  return(plt)


}

plot_random_epochs <- function(epochs,n=100)
{

  if('timepoint' %in% names(epochs))
  {
    epochs <- dplyr::select(epochs,ID,eyetrial,timepoint,value)
    epochs <- tidyr::spread(epochs,timepoint,value)
  }


  epochs <- dplyr::select(epochs,matches('^t[_0-9]'))

  tmp <- dplyr::sample_n(epochs,n)
  matplot(t(as.matrix(tmp)),type='l',lty=1)

}


plot.rois <- function(obj,which='all',crosshairs=T){

  df <- rois2df(obj)

  if(which !='all')
    df <- subset(df,name %in% which)

  p <- ggplot2::ggplot() + ggplot2::geom_polygon(data=df, ggplot2::aes(x=x,y=y,group=name),fill='gray',alpha=0.5)+


    ggplot2::geom_text(data=unique(df[c('xcenter','ycenter','name')]),ggplot2::aes(x=xcenter,y=ycenter,label=name),color='white') +
    ggplot2::xlim(c(0,obj$resolution[1])) + ggplot2::ylim(c(obj$resolution[2],0)) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'black'),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())



  if(crosshairs){
    p <- p +
      ggplot2::geom_hline(yintercept=384,color='red',size=0.3,linetype='dashed') +
      ggplot2::geom_vline(xintercept=512,color='red',size=0.3,linetype='dashed')
  }

  p
  return(p)

}


