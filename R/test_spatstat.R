# library(edfR)
# library(spatstat)
# library(ggplot2)
#
# xrange <- c(0,1024)
# yrange <- c(0,768)
#
#
# ev <- edf.events('~/Dropbox/finalized_experiments/capture_endo/data/105_cap.edf')
#
# el <- ellipse(50,100,centre=c(250,600),phi=125)
# plot(as.mask(el))
#
# circ <- disc(radius=50,centre=c(512,384))
#
# ev$hit <- as.numeric(inside.owin(ev$gavx,ev$gavy,circ))
#
# temp <- subset(ev,hit==1)
#
# qplot(gavx,gavy,data=temp,xlim=c(0,1024),ylim=c(0,768))
#
# screen <- owin(xrange,yrange)
# points <- ppp(ev$gavx,ev$gavy,window=screen)
#
#
# endocoords = radialCoords(512,384,12,240,0);
# exocoords = radialCoords(512, 384,6, 280,30);
#
# endocoords[c(2,4,6,8,10,12),] <- exocoords
#
#
# roipoints <- endocoords
#
# rois <- ppp(roipoints[,1],roipoints[,2],
#             window=screen)
#
# plot(rois)
#
#
#
# #testing rotation
# # rois2 <- rois
# # coords(rois2) <- coords(rotate(rois,25,centre='midpoint'))
# # plot(superimpose(rois,rois2))
#
#
# #create a bunch of ellipses
# allels <- list()
# roidfs <- list()
# angles <- rev(seq(0,2*pi,(2*pi)/nrow(coords(rois)))) #angles for the "flower" pattern. Accomodates arbitrary # of ROIs
#
# for(i in 1:nrow(coords(rois))){
#
# eltmp <- spatstat::ellipse(50,100,centre=coords(rois)[i,],phi=angles[[i]])
#
# allels[[i]] <- eltmp
#
# eltmp <- as.data.frame(vertices(eltmp))
# eltmp$roi <- i
# eltmp$xcenter <- coords(rois)[i,1]
# eltmp$ycenter <- coords(rois)[i,2]
#
# roidfs[[i]]<- eltmp
#
# }
#
# ## combine into 1 mask
# # comb <- as.mask(do.call('union.owin',as.list(allels)),dimyx=c(768,1024))
#
# #or do as a layered object
# # allelipse <- layered(LayerList = allels)
# # plot(allelipse)
#
#
# # test <- layered(rois = comb,points=points,plotargs=list(list(col='gray'),list(col='blue',pch=16,cex=0.5)))
# # iplot(test)
#
#
# xy <- as.data.frame(points)
# test <- do.call('rbind',roidfs)
#
# ## plotting just the ROIs as polygons
# # qplot(x,y,data=test,group=roi,fill=factor(roi),geom='polygon')
#
# ggplot() + geom_polygon(data=test, aes(x=x,y=y,group=roi),fill='gray',alpha=0.5)+
#   geom_point(data=coords(points),aes(x=x,y=y),color='yellow',size=.7) +
#   geom_text(data=test,aes(x=xcenter,y=ycenter,label=roi),color='white') +
#   geom_hline(yintercept=384,color='red',size=0.3,linetype='dashed') +
#   geom_vline(xintercept=512,color='red',size=0.3,linetype='dashed') +
#   xlim(c(0,1024)) + ylim(c(768,0)) +
#   theme(panel.background = element_rect(fill = 'black'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
#
#
#
#
#
# test <- layered(p=points,r=rois)
#
# plot(test)
#
#
# #applying a shift to all fixations.
# testshift <- spatstat::shift(points,c(20,-20))
# Window(testshift) <- screen
# plot(testshift,pch=16)
# plot(points,add=T,pch=16,col='red') #dropshadow effect!
#
# #applying a rotation ot all fixations
# testrotate <- rotate(points,angle = 10,centre='centroid')
# Window(testrotate) <- screen
# plot(testrotate,pch=16)
# plot(points,add=T,pch=16,col='red') #dropshadow effect!
#
#
#
# # plot(subset(points,V5<=100))
#
# #indicating
# marks(points) <- cbind(as.factor(inside.owin(points,w=allels[[5]])),as.factor(inside.owin(points,w=allels[[10]])))
#
# plot(split(points,f='V1'))
