
library(raster)

#get observations
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'

obs<-read.csv(paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear2015-07-12.csv"))

month<-c("01","02","03","04","05","06","07","08","09","10","11","12")

monthdat<-list()
for(i in 1:length(month)){
  monthObs<-subset(obs, month==i)
  #get spatial data
  tmax<-raster(paste0("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTClimate_mmn/tmax/eMAST_ANUClimate_mmn_tmax_v1m0_",
                      month[i],".nc"),band=1)
  #make spatial points
  IDLocs <- SpatialPoints(cbind(monthObs$lon, monthObs$lat),
                          proj4string = CRS('+proj=longlat +datum=WGS84'))
  
  IDclim<-data.frame(extract(tmax,IDLocs))
  dat<-cbind(monthObs,IDclim)
  colnames(dat)[ncol(dat)]<-"mmtmax"
  monthdat[[i]]<-dat
}

findat<-do.call("rbind",monthdat)

write.csv(findat,paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear",as.Date(Sys.Date()),'.csv'),row.names=FALSE)

                        
                        


