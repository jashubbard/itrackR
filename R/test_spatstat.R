library(edfR)
library(spatstat)
library(ggplot2)

xrange <- c(0,1024)
yrange <- c(0,768)


ev <- edf.events('/Volumes/Data/data_archive/switchrate_eug/101_sr.edf')

el <- ellipse(50,100,centre=c(250,600),phi=125)
plot(as.mask(el))

circ <- disc(radius=50,centre=c(512,384))

ev$hit <- as.numeric(inside.owin(ev$gavx,ev$gavy,circ))

temp <- subset(ev,hit==1)

qplot(gavx,gavy,data=temp,xlim=c(0,1024),ylim=c(0,768))

points <- ppp(ev$gavx,ev$gavy,xrange,yrange)

rois <- ppp(c(512,612,686,712,686,612,512,412,339,312,339,412),
            c(584,558,484,384,285,211,184,211,284,384,484,558),
            xrange, yrange)

plot(rois)


#testing rotation
# rois2 <- rois
# coords(rois2) <- coords(rotate(rois,25,centre='midpoint'))
# plot(superimpose(rois,rois2))


#create a bunch of ellipses
allels <- list()
angles <- rev(seq(0,2*pi,(2*pi)/nrow(coords(rois)))) #angles for the "flower" pattern. Accomodates arbitrary # of ROIs

for(i in 1:nrow(coords(rois))){

allels[[i]] <- spatstat::ellipse(50,100,centre=coords(rois)[i,],phi=angles[[i]])

}

#combine into 1 mask
comb <- as.mask(do.call('union.owin',as.list(allels)),dimyx=c(768,1024))

#or do as a layered object
allelipse <- layered(LayerList = allels)
plot(allelipse)


test <- layered(rois = comb,points=points,plotargs=list(list(col='gray'),list(col='blue',pch=16,cex=0.5)))
iplot(test)


xy <- coords(points)

library(grid)

img <- rasterGrob(t(as.matrix(comb)),interpolate=T)
ggplot(xy,aes(x,y)) +

  qplot(1:10, 1:10, geom="blank") + geom_raster(aes(x,y,value),data=test,color='gray') +xlim(c(0,1024)) + ylim(c(0,768))
# marks(points) <- crossdist(points,rois)




test <- layered(p=points,r=rois)

iplot(test)

# plot(subset(points,V5<=100))

iplot(points)
marks(points) <- as.factor(inside.owin(points,w=allels[[10]]))

plot(split(points))
