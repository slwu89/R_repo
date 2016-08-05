#################################################################
######Experiments in d3.js through R for Data Visualization######
######Sean Wu 8/4/2016###########################################
#################################################################

library(threejs)
library(rgdal)

memory.limit(size=540000)

globejs(img="C:/Users/Administrator/Dropbox/GitHub/R_repo/world.jpg",atmosphere=TRUE)

# globejs(img=system.file("images/world.jpg", package="threejs"),atmosphere=TRUE)

#sample with world temperature data
x <- readGDAL("C:/Users/Administrator/Dropbox/GitHub/R_repo/worldTemp.TIFF")
x <- as.data.frame(cbind(coordinates(x), x@data[,1]))
names(x) <- c("long","lat","value")
x <- x[x$value < 255,]
x <- na.exclude(x)

col = topo.colors(length(unique(x$value)))[x$value]


globejs(img="C:/Users/Administrator/Dropbox/GitHub/R_repo/worldTemp.jpg",
        lat=x$lat,long=x$long,val=x$value,color=col,pointsize=0.5,atmosphere=TRUE)


globejs(lat=x$lat,long=x$long,val=x$value,color=col,pointsize=0.5,atmosphere=TRUE)
