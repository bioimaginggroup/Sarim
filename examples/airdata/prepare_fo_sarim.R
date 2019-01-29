data<-data.final.p1
data<-data[data$timestamp>="2019-01-21 00:00:00"]
data$lat<-round(1000*data$lat,0)
data$lon<-round(500*data$lon,0)
latgrid <- seq(min(data$lat),max(data$lat),by=1)
longrid <- seq(min(data$lon),max(data$lon),by=1)
lonlatid<-data.frame()
vecdat<-data.frame()
id<-0
for (i in latgrid)
  for (j in longrid)
  {
    id<-id+1
    lonlatid<-rbind(lonlatid,c(id,i,j))
  }

names(lonlatid)<-c("id","lat","lon")
for (i in 1:dim(data)[1])
{
  d<-data[i,]
  w<-which((lonlatid$lat==d$lat)&(lonlatid$lon==d$lon))
  vecdat<-rbind(vecdat,c(lonlatid$id[w],d$value,d$temp,d$trafficvol,d$humi))
}
names(vecdat)<-c("id","value","temp","trafficvol","humi")
vecdat<-vecdat[!is.na(vecdat$value),]

Q1<-BayGMRF::makeQ(length(latgrid),"rw1")
Q2<-BayGMRF::makeQ(length(longrid),"rw1")
I1<-diag(length(latgrid))
I2<-diag(length(longrid))

Q<-kronecker(Q1,I2)+kronecker(I1,Q2)
save.image("examples/airdata/prep.Rdata")
