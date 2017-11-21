#Calculate the quantiles on circular data
#make graph of Australia wide data

rm(list = ls())
library(raster)
library(circular)
#library(CircStats)
library(gtools)
library(Hmisc)

CalcbinSize <- function(epochDates){
  d <- density(epochDates,na.rm = TRUE,from = 1, to = 366)
  ap <- approxfun(x=d$x, y=d$y)
  binSize<-integrate(ap, 1, 30.5)[[1]]
  for (mnth in 1:11){
    binSize<-c(binSize,integrate(ap, 30.5*mnth, (30.5*mnth)+30.5)[[1]])
  } 
  return(binSize)
}

calcskew <- function(circulardates){
  datesC<-circular( circulardates,modulo ="2pi", units="radians", rotation="counter")
  Rbar<-rho.circular(datesC)#average clustering 0 is uncluseted 1 is all at same location
  V<-1-Rbar#sample circular variance
  t2t<- trigonometric.moment(datesC, p=2, center=TRUE)
  bbar2 <- t2t$sin
  skewness <- bbar2/(V**(3/2)) #skewness 
  return(round(skewness,2))
}

#Koeppen zones
#   41 Equatorial
#   35 Tropical
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

#working directory
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
dat<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2015-10-07.csv'))

#figure out data types
#egg
dat$type<- with(dat,ifelse((startEgg==startEgg & is.na(startHatch) & is.na(startYoung) & is.na(endBuild)),"solitary-egg",type))
#unknown
dat$type<- with(dat,ifelse((is.na(startEgg)&is.na(startEgg) & is.na(startHatch) & is.na(startYoung) & is.na(endBuild)),"unknown",type))
#young
dat$type<- with(dat,ifelse((is.na(startEgg)&is.na(startEgg) & is.na(startHatch) & !is.na(startYoung) & is.na(endBuild)),"young",type))
# multi
dat$type<- with(dat,ifelse((!is.na(startHatch)),"multi",type))     
dat$type<- with(dat,ifelse((!is.na(startEgg) &!is.na(startEgg) & !is.na(startYoung)),"multi",type))     
dat$type<- with(dat,ifelse((!is.na(startEgg) &!is.na(startEgg) & !is.na(endBuild)),"multi",type))
dat$type<- with(dat,ifelse((is.na(startEgg) &is.na(startEgg) & !is.na(endBuild)& !is.na(startYoung)),"multi",type))
#build plus unknown
dat$type<- with(dat,ifelse((is.na(startEgg) & is.na(startEgg) & !is.na(endBuild) & is.na(startYoung) & is.na(startHatch)& !is.na(startUnknown)),"unknown",type)) 


#add koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen, rclmat)


dat$koeppen<-extract(koeppen,data.frame(cbind(dat$lon, dat$lat)))
#convert data to radians to make circular vector
dat$Radians <-(dat$DOY_PL/366*360)*pi / 180
#get list of incubation, fledging, clutch size etc
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#get list of species not to include
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
dat<-dat[dat$Scientific.Name %nin% RMspecies,]
dat<-subset(dat,!is.na(dat$koeppen))
inc<-inc[inc$ScientificName %nin% RMspecies,]

species<-inc$ScientificName

dat<-subset(dat,year!= 1770 | Scientific.Name != "Passer domesticus")
#write.csv(dat,paste0(dat.dir,'PointOfLayDayOfYear2016-09-20.csv'))


#species<-"Pomatostomus ruficeps"
#loop through species and convert to circular and find quantiles
quantdat<-list()
peakBreeding<-list()
for (i in 1:length (species)){
  spdat<-subset(dat,Scientific.Name==species[i])#get species data
  CommonName<-as.character(subset(inc,ScientificName==species[i],"Common.me",drop=TRUE))
  Order<-as.character(subset(inc,ScientificName==species[i],order,drop=TRUE))
  #########################
  # koeppen zone analyses
  #########################  
  #koeppen zones
  koep<-c(35,32,22,13,3)
  kpName<-c('Tropical','Subtropical','Desert','Grassland','Temperate')
  koepObs<-list()
  #colour<-c('darkgreen','green','chartreuse','darkgoldenrod4','purple','blue3')
  koepQuant<-list()
  for(ii in 1:length(koep)){
    kpdat<-as.numeric(subset(spdat,koeppen==koep[ii],Radians)[,1])
    if(length(kpdat)>=50){
      #date of quantiles as long as over 50 obs in a region, number comes from Joys and Cricks
      kpquant<-quantile.circular(kpdat,c(.05,.5,.95),type=8)
      kp5<-round((kpquant[[1]]*180/pi)/360*365) 
      kp50<-round((kpquant[[2]]*180/pi)/360*365)
      kp95<-round((kpquant[[3]]*180/pi)/360*365)
      #Breeding Period
      kpStartDate<-as.numeric(kp5)
      kpEndDate<-as.numeric(kp95)
      kpBPL<-ifelse (kpEndDate < kpStartDate, 365-kpStartDate+kpEndDate,kpEndDate - kpStartDate )
      #get rbar and skewness for all data and 90% of data
      kpdatC<-circular(kpdat,modulo ="2pi", units="radians", rotation="counter") #to circular
      kpRbar<-round(rho.circular(kpdatC),2)#average clustering 0 is uncluseted 1 is all at same location
      kpSkew<-calcskew(kpdat)
      #90 percent of radians
      kpdat90<-if (kpquant[[1]]<kpquant[[3]]){
        c(kpdat[kpdat>kpquant[[1]] & kpdat<kpquant[[3]]])
      } else c(kpdat[kpdat>kpquant[[1]]],kpdat[kpdat<kpquant[[3]]])
      kpdat90C<-circular(kpdat90,modulo ="2pi", units="radians", rotation="counter")
      #plot(kpdat90C)
      kpRbar90<-round(rho.circular(kpdat90C),2)
      kpSkew90<-calcskew(kpdat90)
      #put together
      kpSummary<-cbind(paste0(kpName[ii]),paste(species[i]),CommonName,Order,length(kpdat),kp5,kp50,kp95,kpBPL,kpRbar,kpSkew,kpRbar90,kpSkew90)
      colnames(kpSummary)<-c('Region','Species','CommonName','Order','ObservationCount','Quantile5','Quantile50','Quantile95','BreedingPeriod','RbarAll','SkewAll','Rbar90','Skew90')
      
      ###90% of data, remove obs outside of peak breeding
      kpObs90<-if (kpquant[[1]]<kpquant[[3]]){
      subset(spdat, Radians >=kpquant[[1]] & Radians <= kpquant[[3]] & koeppen==koep[ii])
      }else rbind(subset(spdat, Radians >=kpquant[[1]] &koeppen==koep[ii]),subset(spdat, Radians <=kpquant[[3]]&koeppen==koep[ii]))
      
      koepObs[[ii]]<-kpObs90
      koepQuant[[ii]]<-kpSummary
      rm(kmSummary)
    } 
  }
  quantdat[[i]]<-if (length(koepQuant)>0) do.call("rbind",koepQuant)
  peakBreeding[[i]]<-if (length(koepObs)>0)do.call("smartbind",koepObs)
}
  

alldat<-na.omit(do.call("rbind",quantdat)[,c('Region','Species','CommonName','Order',
                          'ObservationCount','Quantile5','Quantile50','Quantile95','BreedingPeriod','RbarAll','SkewAll','Rbar90','Skew90')])
peakBreedingAll<-do.call("rbind",peakBreeding)

write.csv(peakBreedingAll,paste0(dat.dir,'PeakBreedingPointOfLayDayOfYear',
          as.Date(Sys.Date()),'.csv'),row.names=FALSE)

write.csv(alldat,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles',
                        as.Date(Sys.Date()),'.csv'),row.names=FALSE)















# # =====================================
# 
# # circular plot of Australia with density curves for biomes
# 
# # ===================================== 
# #subset peakbreeding to keep species I want
# 
# dat<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-07-12.csv")
# 
# inc2<-subset(inc,remove!=1)
# 
# species2<-inc2$ScientificName
# 
# PBAll<-peakBreedingAll[peakBreedingAll$Scientific.Name %in% species2,]
# 
# PBalldat<-alldat[alldat$Species %in% species2,]
# 
# 
# PLc<- circular(PBAll$Radians, units = "radians", zero=pi/2, 
#              rotation = "clock")
# 
# #get midpoints of 12 months on circle
# bins<-12
# arc <- (2 * pi)/bins
# sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
# midpoints <- circular(seq(arc/2, 2 * pi - pi/bins, length = bins), units = "radians", zero=pi/2, 
#                       rotation = "clock")
# 
# #plot(PLc, cex=1.1, bin=720, stack=TRUE, sep=0.035, shrink=1.8)
# rose.diag(PLc, bins=12, col="darkgrey", cex=1.1, prop=1.7,axes = FALSE,ticks = FALSE, shrink=1.5)
# #rose.diag(PLc, bins=12, col="darkgrey", cex=1.1, prop=1.3,rotation='clock',zero=pi/2,axes = FALSE,ticks = FALSE, shrink=1.8)
# axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)
# #ticks.circular(midpoints, rotation='clock', tcl=0.075)
# #lines(density.circular(PLc, bw=40), lwd=2, lty=1)#Australia Density
# 
# koepClass<-c(35,32,22,13,3)
# 
# #Koeppen zones
# #   41 Equatorial
# #   35 Tropical
# #   32 Subtropical
# #   22 Desert
# #   13 Grassland
# #   3 Temperate
# 
# kpColour<-c("chartreuse4","darkmagenta", "darkred","darkgoldenrod1","cornflowerblue")
# 
# for (k in 1:length(koepClass)){
#   kpdat<-circular(subset(dat,koeppen==koepClass[k],select = Radians), units = "radians", zero=pi/2, 
#                   rotation = "clock")
#   lines(density.circular(kpdat, bw=40), lwd=3, lty=k,col = kpColour[k])
#   arrows.circular(mean(kpdat), lwd=3, lty=k,col = kpColour[k])
# }
# 
# legend("topright",legend=c("Desert","Grassland","Subtropical","Temperate","Tropical"),
#        pch=1 , lwd=3, bty="n", text.font=3, 
#        col = c("darkred","darkgoldenrod1","darkmagenta","cornflowerblue","chartreuse4"))
#          
#          

# =====================================

# rose diagram of the number of birds species breeding and 
# kernal density lines of number of bird species breeding in each biome

# ===================================== 

#get unique species per day for Australia

#species 1 obs per a day
#dat$MOY_PL<-ceiling(dat$DOY_PL/7.038462)
PBAll$MOY_PL<-PBAll$DOY_PL
Ausdat<-PBAll[!duplicated(PBAll[c("Scientific.Name","MOY_PL")]),]
PLc<- circular(Ausdat$Radians, units = "radians", zero=pi/2, 
               rotation = "clock")

#get midpoints of 12 months on circle
bins<-52
arc <- (2 * pi)/bins
sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
midpoints <- circular(seq(arc/2, 2 * pi - pi/bins, length = bins), units = "radians", zero=pi/2, 
                      rotation = "clock")

monthlymidpoints <- circular(seq(((2 * pi)/13)/2, 2 * pi - pi/13, length = 13), units = "radians", zero=pi/2, 
                      rotation = "clock")
rose.diag(PLc, bins=52, col="darkgrey", cex=1.1, prop=4.5,axes = FALSE,ticks = FALSE, shrink=1.5)
wk<-round((seq(((2 * pi)/13)/2, 2 * pi - pi/13, length = 13)*180/pi)/360*52)
axis.circular(at=monthlymidpoints, labels=wk, cex=1.1,tcl=0.075)
koepClass<-c(35,32,22,13,3)

kpColour<-c("chartreuse4","darkmagenta", "darkred","darkgoldenrod1","cornflowerblue")

for (k in 1:length(koepClass)){
  kpdat<-subset(PBAll,koeppen==koepClass[k])
  kpdat<- kpdat[!duplicated(kpdat[c("Scientific.Name","MOY_PL")]),]
  kpdat<-circular(kpdat$Radians, units = "radians", zero=pi/2, 
                  rotation = "clock")
  
  
  lines(density.circular(kpdat, bw=10), lwd=3, lty=k,col = kpColour[k])
}
legend("topright",legend=c("Desert","Grassland","Subtropical","Temperate","Tropical"),
       pch=1 , lwd=3, bty="n", text.font=3, 
       col = c("darkred","darkgoldenrod1","darkmagenta","cornflowerblue","chartreuse4"))



# ===================================== 

# Table of Biome Summary

# ===================================== 

head(PBalldat)
regions<-unique(PBalldat$Region)
BSummary<-as.data.frame(paste(c("Australia", kpName)))
colnames(BSummary)<-"Region"

RegionSummary<-list()
for (r in 1:length(regions)){
  D5th<-subset(PBalldat,Region==regions[r],select=Quantile5)
  D95th<-subset(PBalldat,Region==regions[r],select=Quantile95)
  D50th<-subset(PBalldat,Region==regions[r],select=Quantile50)
  vals<-subset(PBalldat,Region==regions[r],select=BreedingPeriod)
  birds<-length(vals[,1])#number of breedin birds
  meanPeriod<-mean(as.numeric(vals[,1]))#average length of breeding period
  sdPeriod<-sd(as.numeric(vals[,1])) #sd of breeding perid
  D50thr <-(as.numeric(na.omit(D50th)[,1])/365*360)*pi / 180#median date to radians
  D50thrC<-circular(D50thr,modulo ="2pi", units="radians", rotation="counter") #to circular
  #plot(D50thrC)
  AvgMedR<-mean(circular(D50thr,modulo ="2pi", units="radians", rotation="counter"))
  AvgMed<-round((AvgMedR[[1]]*180/pi)/360*365)#average median date as day
  Rbar<-rho.circular(D50thrC)#average clustering 0 is uncluseted 1 is all at same location
  #message(Rbar)
  V<-1-Rbar#sample circular variance
  #change this to average skewness across all species
  #skewness, need at least 3 data points to calc
  #negative numbers skewed counter clockwise, poitive clockwise, 0 = not skewed
  t2t <- trigonometric.moment(D50thrC, p=2, center=TRUE)
  bbar2 <- t2t$sin
  skewness <- bbar2/(V**(3/2)) #skewness 
  #abar2 <- t2t$cos
  #hatk <- (abar2-Rbar**4)/(V**2) ; #kurtosis
  region<-regions[r]
  RegionSummary[[r]]<-cbind(region,birds,meanPeriod,sdPeriod,AvgMed,skewness,Rbar)
}

RegionSummaryFinal<-as.data.frame(do.call("rbind",RegionSummary))#convert to data.frame
write.csv(RegionSummaryFinal,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/RegionalSummary',as.Date(Sys.Date()),'.csv'),row.names=FALSE)





