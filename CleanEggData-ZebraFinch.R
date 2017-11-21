


#clean and combine breeding observation data from all museums, ATLAS, Nest Record Scheme
#seperate epoch dat into nest building, egg, young, and unknown event
#fix lat. values that are not negative
#check observations falll within Australia based on eMAST data
#calculate epoch date since 1970-1-1
#write table of species with more than 100 breeding observations


rm(list = ls())

library(raster)
library(car) 
library(stringr)
library(Hmisc)
library(lubridate)
library(gtools)

indir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/egg'

###########################################################################
######NT data
NT<-read.csv(paste0(indir,'/NT_Museum/MAGNT_bird eggs_Jan2015.csv'))
NT<-subset(NT, !is.na(Lat_degrees) & !is.na(Long_degrees))
NT$Lat_minutes <- replace(NT$Lat_minutes, is.na(NT$Lat_minutes),0)
NT$Long_minutes <- replace(NT$Long_minutes, is.na(NT$Long_minutes),0)
NT$lat<-NT$Lat_degrees + (NT$Lat_minutes/60)
NT$long<-NT$Long_degrees + (NT$Long_minutes/60)
NT$Sci.Name<-paste(NT$Genus,NT$Species,sep=' ')
NT$Source <- 'NorthernTerritoryMuseumAndArtGallery'
NT$Coll_Month<-as.numeric(recode(NT$Coll_Month,"'Jan'=1;'Feb'=2;'Mar'=3;'Apr'=4;
                                 'May'=5;'Jun'=6;'Jul'=7;'Aug'=8;'Sep'=9;'Oct'=10;
                                 'Nov'=11;'Dec'=12")) 
NT$epoch<-as.numeric(as.Date(paste0(NT$Coll_Year,'-',NT$Coll_Month,'-',NT$Coll_Day)))
NT<-subset(NT, !is.na(epoch))
NT$startegg<-NT$epoch
NT$endegg<-NT$epoch
NT$starthatch<-NA
NT$startyoung<-NA
NT$endbuild<-NA
NT$startUnknown<-NA
NT$endUnknown<-NA
NT$type<-'solitary-egg'
eggdat<-NT[c('Sci.Name','lat','long', 'Source','startegg',
             'endegg','type')]#,'starthatch','startyoung','endbuild','startUnknown','endUnknown')]
colnames(eggdat)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg',
                     'endEgg','type')#,'startHatch','startYoung','endBuild','startUnknown','endUnknown')

###########################################################################
######South Australia Museaum data
SA<-read.csv(paste0(indir,'/SA_Museam/SA_egg_records.csv'))
SA$lat_min <- replace(SA$lat_min, is.na(SA$lat_min),0)
SA$long_min <- replace(SA$long_min, is.na(SA$long_min),0)
SA<-subset(SA, !is.na(lat_hour) & !is.na(lat_hour))
SA$lat<-SA$lat_hour + (SA$lat_min/60) 
SA$long<-SA$long_hour + (SA$long_min/60)
SA$Sci.Name<-SA$Taxon
SA$Source <- 'SouthAustraliaMuseum'
SA$epoch<-as.numeric(as.Date(paste0(SA$Year,'-',SA$Month,'-',SA$Day)))
SA<-subset(SA, !is.na(epoch))
SA$startegg<-SA$epoch
SA$endegg<-SA$epoch
SA$type<-'solitary-egg'
SA<-SA[c('Sci.Name','lat','long', 'Source','startegg','endegg','type')]
colnames(SA)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')

eggdat<-smartbind(eggdat,SA)

###########################################################################
######Queensland
QU<-read.csv(paste0(indir,'/QL_Museum/egg_nest_sepLATLONG.csv'))
QU<-subset(QU,Field.Coll.Specimen.Category=="Egg(s)" |
             Field.Coll.Specimen.Category=="Egg(s), Nest" |
             Field.Coll.Specimen.Category== "Egg(s), Skeletal parts Nest" |
             Field.Coll.Specimen.Category== "Spirit, Egg(s)")
QU$lat_min <- replace(QU$lat_min, is.na(QU$lat_min),0)
QU$long_min <- replace(QU$long_min, is.na(QU$long_min),0)
QU$lat_sec <- replace(QU$lat_sec, is.na(QU$lat_sec),0)
QU$Long_sec <- replace(QU$Long_sec, is.na(QU$Long_sec),0)
QU$lat<-QU$lat_deg + (QU$lat_min/60) + (QU$lat_sec/3600) *-1
QU$long<-QU$long_deg + (QU$long_min/60)+ (QU$Long_sec/3600) *-1
QU$Sci.Name<-paste(QU$Genus, QU$Species, sep=" ")
QU$Source <- 'QueenslandMuseum'
QU$type<-'solitary-egg'
QU$epoch<-as.numeric(as.Date(paste0(QU$year,'-',QU$month,'-',QU$day)))
QU<-subset(QU, !is.na(epoch))
QU$startegg<-QU$epoch
QU$endegg<-QU$epoch
QU<-QU[c('Sci.Name','lat','long', 'Source','startegg','endegg','type')]
colnames(QU)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')

eggdat<-smartbind(eggdat,QU)


###########################################################################
######Australia Museaum data
ausmus<-read.csv(paste0(indir,'/Australia Museum/AUSMUS_egg_records.csv'))
ausmus<-subset(ausmus, !is.na(start_latitude) & !is.na(start_longitude) & !is.na(year) 
               & !is.na(month)& !is.na(day))
ausmus$Sci.Name<-paste(ausmus$genus,ausmus$species,sep=' ')
ausmus$Source <- 'AustraliaMuseum'
ausmus$epoch<-as.numeric(as.Date(paste0(ausmus$year,'-',ausmus$month,'-',ausmus$day)))
ausmus<-subset(ausmus, !is.na(epoch))
ausmus$startegg<-ausmus$epoch
ausmus$endegg<-ausmus$epoch
ausmus$type<-'solitary-egg'
ausmus$precision<-ausmus$coordinateprecision
ausmus$precision<-recode(ausmus$precision,"'1km-10km'='10'; '10km-100km'='100'; '10m-100m' = '.1';'100m-1km'='1'")
ausmus<-ausmus[c('Sci.Name','start_latitude','start_longitude','Source','startegg','endegg','type','precision')]
colnames(ausmus)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type','precision')
eggdat<-smartbind(eggdat,ausmus)

###########################################################################
######Victoria

Vic1<-read.csv(paste0(indir,'/Vic_Museum/OZCAM/data.csv'))
Vic1$eggcode <- lapply(list(c(1,2),c(3,4)),function(i) substr(Vic1$Catalog.Number,i[1],i[2]))[[1]]
#select out egg obs
Vic1<-subset(Vic1,eggcode=='BE')
#Vic1<- subset(Vic1, Country...parsed == 'Australia')
Vic1<-Vic1[,c('Scientific.Name','Latitude...processed','Longitude...processed','Year...parsed','Month...parsed','day')]
colnames(Vic1)<-c('Scientific.Name','lat', 'lon','year','month','day')
Vic1$sourceName<-'MuseumVictoria'

Vic2<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection1.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..','year','month','Day')]
colnames(Vic2)[7]<-'day'
Vic3<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection2.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..','year','month','day')]
Vic4<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection3.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..','year','month','day')]
Vic5<-read.csv(paste0(indir,'/Vic_Museum/MVEgg Collection4.csv'))[c('Genus...Version.1.2.elements..1..', 'Species...Version.1.2.elements..2..','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..','year','month','day')]
NVic<-rbind(Vic2,Vic3,Vic4,Vic5)
NVic$Species<-paste(NVic$Genus...Version.1.2.elements..1..,NVic$'Species...Version.1.2.elements..2..')
NVic$sourceName<-'MuseumVictoria'
NVic<-subset(NVic, !is.na('Latitude...Version.1.2.elements..3..') 
             & !is.na(Longitude...Version.1.2.elements..3..) 
             & !is.na(year) & !is.na(month)& !is.na(day))

NVic$epoch<-as.numeric(as.Date(paste0(NVic$year,'-',NVic$month,'-',NVic$day)))
NVic<-subset(NVic, !is.na(epoch))
NVic$startegg<-NVic$epoch
NVic$endegg<-NVic$epoch
NVic$type<-'solitary-egg'
NVic<-NVic[c('Species','Latitude...Version.1.2.elements..3..','Longitude...Version.1.2.elements..3..'
             ,'sourceName','startegg','endegg','type')]
colnames(NVic)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')
eggdat<-smartbind(eggdat,NVic)

###########################################################################
######ANWC
ANWC<-read.csv(paste0(indir,'/ANWC_OZCAM/OZCAM_eggs.csv'))
ANWC<-subset(ANWC,Country...parsed=='Australia')
ANWC$Source<-'AustralianNationalWildlifeCollection'
ANWC<-subset(ANWC, !is.na(Latitude...processed) 
             & !is.na(Longitude...processed) 
             & !is.na(Year...parsed) & !is.na(Month...parsed)& !is.na(day))
ANWC$epoch<-as.numeric(as.Date(paste0(ANWC$Year...parsed,'-',ANWC$Month...parsed,'-',ANWC$day)))
ANWC<-subset(ANWC, !is.na(epoch))
ANWC$startegg<-ANWC$epoch
ANWC$endegg<-ANWC$epoch
ANWC$type<-'solitary-egg'
ANWC<-ANWC[c('Scientific.Name','Latitude...processed','Longitude...processed',
             'Source','startegg','endegg','type')]
colnames(ANWC)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')
eggdat<-smartbind(eggdat,ANWC)

###########################################################################
###### Western Australia
WA<-read.csv(paste0(indir,'/WA_Museum/WA_museum.csv'))
WA$precision<-ifelse(!is.na(WA$Coordinate.Uncertainty.in.Metres...parsed),
                     WA$Coordinate.Uncertainty.in.Metres...parsed/1000,"NA")
#get the ID for egg records
rego<-read.csv(paste0(indir,'/WA_Museum/eggs-regno.csv'))[,1]
WAeggs<-WA[WA$'Catalog.Number' %in% rego,] 
WAeggs<-subset(WAeggs, Country...parsed=='Australia')
WAeggs$Source<-'WesternAustraliaMuseum'
WAeggs<-subset(WAeggs, !is.na(Latitude...processed) 
               & !is.na(Longitude...processed) 
               & !is.na(Year...parsed) & !is.na(Month...parsed)& !is.na(day))
WAeggs$epoch<-as.numeric(as.Date(paste0(WAeggs$Year...parsed,'-',
                                        WAeggs$Month...parsed,'-',WAeggs$day)))
WAeggs<-subset(WAeggs, !is.na(epoch))
WAeggs$startegg<-WAeggs$epoch
WAeggs$endegg<-WAeggs$epoch
WAeggs$type<-'solitary-egg'
WAeggs<-WAeggs[c('Scientific.Name','Latitude...processed','Longitude...processed',
                 'Source','startegg','endegg','type','precision')]
colnames(WAeggs)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type','precision')
eggdat<-smartbind(eggdat,WAeggs)


###########################################################################
##### Nest Record Scheme
NRS<-read.csv(paste0(indir,'/BLA_NRS/NRSExtract.csv'))
NRS$YEAR<- ifelse(NRS$YEAR>2015,NA,NRS$YEAR)
NRS$MONTH<- ifelse(NRS$MONTH>12,NA,NRS$MONTH)
NRS$DAY<- ifelse(NRS$DAY>31,NA,NRS$DAY)
NRS<-with(NRS,NRS[!is.na(DAY) & !is.na(MONTH)  & !is.na(YEAR) &
                    !is.na(date) & !is.na(Lat) & !is.na(Lon),])
NRS$ID<-1:nrow(NRS)

#code for E=eggs, H = eggs hatching,B = nest being built,U = Unknown
E<-with(NRS,NRS[NRS_VISt_Event=='M' 
                | NRS_VISt_Event == 'm' 
                | lkpEvent_Event =='changeover/both off'
                | lkpEvent_Event =='male at nest'
                | NRS_VISt_Event=='F' 
                | NRS_VISt_Event == 'f' 
                | lkpEvent_Event =='female at nest'
                | lkpEvent_Event =='cold eggs'
                | NRS_VISt_Event=='u' 
                | NRS_VISt_Event == 'U'
                | NRS_VISt_Event=='x' 
                | NRS_VISt_Event == 'X'
                | lkpEvent_Event =='sex unknown at nest'
                | lkpEvent_Event =='warm eggs',])
H<-with(NRS,NRS[NRS_VISt_Event=='H' 
                | lkpEvent_Event == 'eggs hatching',])
B<-with(NRS,NRS[NRS_VISt_Event=='b' 
                |NRS_VISt_Event=='B' 
                |lkpEvent_Event == 'nest being built',])
U<-NRS[NRS$ID %nin% (c(B$ID,H$ID,E$ID)),]
E$EventCode<-'E'
H$EventCode<-'H'
B$EventCode<-'B'
U$EventCode<-'U'
NRS<-rbind(E,H,B,U) 

#if eggs are in nest give EventCode E as long as not already H
df<-subset(NRS, EGGS != '')
df<-with(df,df[EGGS !=0
               & EGGS != '0?'
               & EGGS != 'B'
               & EGGS != 'O'
               & EventCode != 'H',])
fix<-NRS[NRS$ID %in% (df$ID),]
fix$EventCode<-'E'
nofix<-NRS[NRS$ID %nin% (df$ID),]
NRS2<-rbind(fix, nofix)

#if Young are in nest give EventCode Y as long as not already H
df<-subset(NRS2, Y_IN != '')
df<-with(df,df[Y_IN !=0
               & Y_IN != '0?'
               & Y_IN != '?'
               & Y_IN != 'a'
               & Y_IN != 'O'
               & Y_IN != 'U'
               & Y_IN != 'W'
               & EventCode != 'H',])
fix<-NRS2[NRS2$ID %in% (df$ID),]
fix$EventCode<-'Y'
nofix<-NRS2[NRS2$ID %nin% (df$ID),]
NRS3<-rbind(fix, nofix)

#if young outside nest add Y_Out 
df<-subset(NRS3, Y_OUT != '')
df<-with(df,df[Y_OUT != 0
               & Y_OUT != '0?'
               & Y_OUT != '?'
               & Y_OUT != 'O'
               & Y_OUT != 'U'
               & Y_OUT != 'W'
               & EventCode != 'H'
               & EventCode != 'Y',])

#if young are inside the nest
fix<-NRS3[NRS3$ID %in% (df$ID),]
fix$EventCode<-'Y'
fix$EventCode2<-'Y_OUT'
nofix<-NRS3[NRS3$ID %nin% (df$ID),]
NRS4<-smartbind(fix, nofix)

#make isodate, number of dates since 1970
NRS4$epochdate<-as.numeric(as.Date(dmy_hms(NRS4$date)))
NRS4<-subset(NRS4, !is.na(epochdate))

#list of the Nest ID
NRS_Id<-unique(NRS4$SPEC_REF)

#loop through each visit and get details from visits
NRSnestObs <- list()
for(i in 1:length(NRS_Id)){
  subNRS<-subset(NRS4,SPEC_REF==NRS_Id[i])
  SUMsubNRS<-subNRS[1,c('SPEC_REF','Scientific_name','Lat','Lon')]
  SUMsubNRS$SPEC_REF<-NRS_Id[i]
  SUMsubNRS$obsCount<-nrow(subNRS)
  SUMsubNRS$startegg<- ifelse ('E' %in% unique(subNRS$EventCode),
                               min(subset(subNRS,EventCode=='E')$epochdate),NA)
  SUMsubNRS$endegg<-ifelse ('E' %in% unique(subNRS$EventCode),
                            max(subset(subNRS,EventCode=='E')$epochdate),NA)
  SUMsubNRS$starthatch<-ifelse ('H' %in% unique(subNRS$EventCode),
                                min(subset(subNRS,EventCode=='H')$epochdate),NA)
  SUMsubNRS$startyoung<-ifelse ('Y' %in% unique(subNRS$EventCode),
                                min(subset(subNRS,EventCode=='Y')$epochdate),NA)
  SUMsubNRS$endyoung<-ifelse ('Y' %in% unique(subNRS$EventCode),
                              pmax(subset(subNRS,EventCode=='Y' & is.na(EventCode2))$epochdate),NA)
  SUMsubNRS$youngOut<-ifelse ('Y_OUT' %in% unique(subNRS$EventCode2),
                              min(subset(subNRS,EventCode2=='Y_OUT')$epochdate),NA)
  SUMsubNRS$endbuild<-ifelse ('B' %in% unique(subNRS$EventCode),
                              max(subset(subNRS,EventCode=='B')$epochdate),NA)
  SUMsubNRS$startUnknown<-ifelse ('U' %in% unique(subNRS$EventCode),
                                  min(subset(subNRS,EventCode=='U')$epochdate),NA)
  SUMsubNRS$endUnknown<- ifelse ('U' %in% unique(subNRS$EventCode),
                                 max(subset(subNRS,EventCode=='U')$epochdate),NA)
  #   ifelse ('Y_OUT' %in% unique(subNRS$EventCode2),
  #           SUMsubNRS$youngOut<-min(subset(subNRS,EventCode2=='Y_OUT')$epochdate),SUMsubNRS$youngOut<-NA)
  NRSnestObs[[i]]<-SUMsubNRS
  message(i)
}

NRS_eggobs<-do.call('rbind',NRSnestObs)
NRS_eggobs$sourceName<-'NestRecordScheme'
#add in observation type
NRS_eggobs$type<-with(NRS_eggobs, ifelse (is.na(startegg) & is.na(starthatch)&is.na(startyoung),
                                          'unknown',NA))
NRS_eggobs$type<-with(NRS_eggobs, ifelse (is.na(startegg) & is.na(starthatch) & is.na(type) & !is.na(startyoung),
                                          'youngOut',type))
NRS_eggobs$type<-with(NRS_eggobs, ifelse (is.na(starthatch) & is.na(startyoung) & startegg == endegg & is.na (type),
                                          'solitary-egg',type))
NRS_eggobs$type<-with(NRS_eggobs, ifelse ( is.na (type),
                                           'multi-egg',type))
NRS_eggobs<-NRS_eggobs[c('SPEC_REF', 'Scientific_name','Lat','Lon','sourceName','startegg',
                         'endegg','starthatch','startyoung','youngOut','endbuild','startUnknown','endUnknown','type')]
colnames(NRS_eggobs)<-c( 'SPEC_REF','Scientific.Name','lat','lon','sourceName','startEgg',
                         'endEgg','startHatch','startYoung','youngOut','endBuild','startUnknown','endUnknown','type')

eggdat<-smartbind(eggdat,NRS_eggobs)
###########################################################################
######### Queen Victoria Museum 
QV<-read.csv(paste0(indir,'/QV_Museum/QVM.csv'))
#QV<-subset(QV,Country...parsed=='Australia')
QV$Source<-'QueenVictoria'
QV<-QV[,c('Scientific.Name','Latitude...processed','Longitude...processed', 'Year...parsed','Month...parsed','day','Source')]
colnames(QV)<-c( 'Scientific.Name','lat','lon','year','month','day','sourceName')

QVegg<-read.csv(paste0(indir,'/QV_Museum/QVM Egg Data.csv'))
QVegg$Source<-'QueenVictoria'
QVegg$Scientific.Name<-paste(QVegg$Genus, QVegg$Species,sep = ' ')
QVegg<-QVegg[,c('Scientific.Name','lat','long', 'year','month','day','Source')]
colnames(QVegg)<-c( 'Scientific.Name','lat','lon','year','month','day','sourceName')

# QVnest<-read.csv(paste0(indir,'/QV_Museum/QVM NEST Data.csv'))
# QVnest$Source<-'QueenVictoria'
# QVnest$Scientific.Name<-paste(QVnest$Genus, QVnest$Species,sep = ' ')
# QVnest<-QVnest[,c('Scientific.Name','lat','long', 'year','month','day','Source')]
# colnames(QVnest)<-c( 'Scientific.Name','lat','lon','year','month','day','sourceName')
QV<-rbind(QV,QVegg)

QV<-subset(QV, !is.na(lat) & !is.na(lon) & !is.na(year) & !is.na(month)& !is.na(day))
QV$epoch<-as.numeric(as.Date(paste0(QV$year,'-',QV$month,'-',QV$day)))
QV<-subset(QV, !is.na(epoch))
QV$startegg<-QV$epoch
QV$endegg<-QV$epoch
QV$type<-'solitary-egg'
QV<-QV[c('Scientific.Name','lat','lon','sourceName','startegg','endegg','type')]
colnames(QV)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')
eggdat<-smartbind(eggdat,QV)

##############################################################################
#TMAG
tmag<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/egg/TMAG_Museum/TMAGeggsnests.csv')
tmag<-subset(tmag, !is.na(Lat) & !is.na(Long) & !is.na(year) & !is.na(month)& !is.na(day))
tmag$epoch<-as.numeric(as.Date(paste0(tmag$year,'-',tmag$month,'-',tmag$day)))
tmag<-subset(tmag, !is.na(epoch))
tmag$startegg<-tmag$epoch
tmag$endegg<-tmag$epoch
tmag$Source<-'TasmanianMuseumArtGallery'
tmag$type<-'solitary-egg'
tmag<-tmag[c('Species.Name','Lat','Long','Source','startegg','endegg','type')]
colnames(tmag)<-c( 'Scientific.Name','lat','lon','sourceName','startEgg','endEgg','type')
eggdat<-smartbind(eggdat,tmag)

##############################################################################
#ATLAS
ATLAS<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/observation/ATLAS/AtlasSurveysObs.csv')
ATLAS<-subset(ATLAS,Breeding==1)
ATLAS$Source2<-paste0('ATLAS',', ',ATLAS$Source)
ATLAS$sameDATE<-(as.numeric(as.Date(dmy_hms(ATLAS$StartDate))))-(as.numeric(as.Date(dmy_hms(ATLAS$FinishDate))))
ATLAS<-subset(ATLAS,sameDATE > -10 & sameDATE<10) #(accuracy of +- 5 days)
#add Scientific names
birdNames<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/BLA_Working_List_v1.1.csv')
ATLAS<-merge(ATLAS,birdNames,by.x = 'AtlasNo',by.y='SpNo',all.x=TRUE)
ATLAS<-subset(ATLAS,!is.na(Taxon.scientific.name))
ATLAS$epochdate<-floor(apply(cbind((as.numeric(as.Date(dmy_hms(ATLAS$StartDate)))),(as.numeric(as.Date(dmy_hms(ATLAS$FinishDate))))), 1,mean, trim = 0))
ATLAS<-subset(ATLAS, !is.na(epochdate))
ATLAS$startUnknown<-ATLAS$epochdate
ATLAS$endUnknown<-ATLAS$epochdate
ATLAS$type<-'unknown'
ATLAS$precision<-recode(ATLAS$Accuracy,"'0'=''; '1'='.1'; '2' = '.5';'3'='5';'4'='50';'9'='';'10'='';'11'=''")
ATLAS<-ATLAS[,c('Taxon.scientific.name','Lat','Lon','Source2','startUnknown','endUnknown','type','precision')]
colnames(ATLAS)<-c( 'Scientific.Name','lat','lon','sourceName','startUnknown','endUnknown','type','precision')
eggdat<-smartbind(eggdat,ATLAS)

##################################################################
#########eBird
ebird<-read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/observation/eBird/ebd_AU_prv_relFeb-2015/ebd_AU_prv_relFeb-2015.txt',
                header=TRUE,sep="\t",strip.white=TRUE)
ebird$keep<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE =='NY'|BREEDING.BIRD.ATLAS.CODE =='NE'|
                                 BREEDING.BIRD.ATLAS.CODE =='ON'|BREEDING.BIRD.ATLAS.CODE =='PE'|
                                 BREEDING.BIRD.ATLAS.CODE =='FL'|BREEDING.BIRD.ATLAS.CODE =='CS' |
                                 BREEDING.BIRD.ATLAS.CODE =='FY'|BREEDING.BIRD.ATLAS.CODE =='N',1,0))
ebird<-subset(ebird,keep==1)
ebird$lat<-ebird$LATITUDE
ebird$lon<-ebird$LONGITUDE
ebird$Scientific.Name<-ebird$SCIENTIFIC.NAME
ebird$sourceName <- 'eBIRD'
ebird$epoch<-as.numeric(as.Date(ebird$OBSERVATION.DATE))


#young in nest
ebird$type<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE == 'NY'|BREEDING.BIRD.ATLAS.CODE == 'CS',
                               'startYoung',NA))
ebird$startYoung<-with (ebird, ifelse (type=='startYoung', ebird$epoch,NA))

#eggs in nest
ebird$type<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE =='NE'|BREEDING.BIRD.ATLAS.CODE =='ON'|
                                 BREEDING.BIRD.ATLAS.CODE =='PE','solitary-egg',type))
ebird$startEgg<-with (ebird, ifelse (type=='solitary-egg', ebird$epoch,NA))
ebird$endEgg<-with (ebird, ifelse (type=='solitary-egg', ebird$epoch,NA))

#Recently Fledged young
ebird$type<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE == 'FL'|BREEDING.BIRD.ATLAS.CODE == 'FY',
                               'youngOut',type))
ebird$youngOut<-with (ebird, ifelse (type=='youngOut', ebird$epoch,NA)) 

#Nest Probable Visiting probable Nest site (primarily hole nesters).
ebird$type<-with(ebird,ifelse (BREEDING.BIRD.ATLAS.CODE == 'N','unknown',type))
ebird$startUnknown<-with (ebird, ifelse (type=='unknown', ebird$epoch,NA))  

ebird<-ebird[c('Scientific.Name','lat','lon', 'sourceName',"type","startYoung","startEgg",
               "endEgg","youngOut","startUnknown")]
eggdat<-smartbind(eggdat,ebird)


##################################################################
#########ABBBS - keep only one unique lat, long, day


abbbsA<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartA_dd.csv"))
abbbsB<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartB_dd.csv"))
abbbs<-rbind(abbbsA,abbbsB)
abbbs<-abbbs[!duplicated(abbbs[,c("SCIENTIFIC_NAME","Day","Month","Year","LAT","LON","BANDER")]),]
abbbs$sourceName<-'ABBBS'
abbbs$epoch<-as.numeric(as.Date(paste0(abbbs$Year,'-',abbbs$Month,'-',abbbs$Day)))
abbbs$type<-with(abbbs, ifelse (AGE == "NESTLING","youngNest",'unknown'))
abbbs$startUnknown<-with (abbbs, ifelse (type=='unknown', abbbs$epoch,NA))  
abbbs$startYoung<-with (abbbs, ifelse (type=='startYoung', abbbs$epoch,NA))
abbbs$lat<-abbbs$LAT
abbbs$lon<-abbbs$LON
abbbs$Scientific.Name<-abbbs$SCIENTIFIC_NAME 
abbbs$sourceName <- 'ABBBS'

abbbs<-abbbs[,c('Scientific.Name','lat','lon', 'sourceName', 'type' ,'startYoung','startUnknown')]
eggdat<-smartbind(eggdat,abbbs)


#################################################################
# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#finish cleaning egdat
eggdat$Scientific.Name<- capitalize(trim(gsub('  ', ' ', eggdat$Scientific.Name,
                                              fixed=TRUE)))#check for extra spaces and capatlize
eggdat<-subset(eggdat, !is.na(lat) & !is.na(lon) & !is.na(Scientific.Name))#make sure lat and longs are given

#make sure all latitudes are negative
latfix<-subset(eggdat, lat > 0)
latfix$lat<-latfix$lat*-1
latgood<-subset(eggdat, lat <= 0)
eggdat<-rbind(latfix, latgood)

#remove duplicates
eggdat<-eggdat[!duplicated(eggdat),]

#make sure years are as expected
syn<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved.csv'), use ==1)
bad<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved.csv'), use ==0)
eggdat<-eggdat[eggdat$Scientific.Name %nin% bad$scn1,]
fix<-eggdat[eggdat$Scientific.Name %in% syn$scn1,]
nofix<-eggdat[eggdat$Scientific.Name %nin% syn$scn1,]#[c(1:8)]
fix<-merge(fix,syn,by.x = 'Scientific.Name',by.y='scn1')
fix$Scientific.Name <-fix$garnett_name
fix<-fix[,c( 'Scientific.Name','lat','lon','sourceName','startEgg',
             'endEgg','startHatch','startYoung','endBuild','startUnknown','endUnknown','SPEC_REF',
             'type','youngOut','precision')]
eggdat<-smartbind(fix,nofix)
eggdat<-eggdat[!duplicated(eggdat),]

#make sure data is within continental Australia
aus<-raster('/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTBiovars/bio_1.asc',
            crs = '+proj=longlat +datum=WGS84')
xy<-cbind(eggdat$lon,eggdat$lat)
eggdat$outsideAustralia<-is.na(extract(aus,xy))
eggdat<-subset(eggdat,outsideAustralia==FALSE)[,c( 'Scientific.Name','lat','lon','sourceName','startEgg',
                                                   'endEgg','startHatch','startYoung','endBuild','startUnknown','endUnknown','SPEC_REF',
                                                   'type','youngOut','precision')]
#put together as much information about species with >=100 obs and HANZAB books
a<-data.frame(table(eggdat$Scientific.Name))
a<-subset(a,Freq>=100)
a<-merge(a,birdNames,by.x='Var1',by.y='Taxon.scientific.name')[c('Var1','Freq','SpNo','Taxon.name',
                                                                 'Family.name','Family.scientific.name','Order','IUCN.Red.List.category',
                                                                 'Scientific.name..C.B.2008.','Scientific.name..IOC.','Scientific.name..Clements.')]
pgDat<-read.csv('/Users/daisy/Google Drive/PhD/birdTraits/Garnett paper/NSD-Data Descriptor/Australian Bird Data Version 1_0 for nest behaviours test_dd.csv')[c("book","pgNu","Taxon.scientific.name..2.")]
a<-merge(a,pgDat,by.x='Var1',by.y='Taxon.scientific.name..2.',all.x=TRUE )
HANZAB<-read.csv('/Users/daisy/Google Drive/PhD/birdTraits/Garnett paper/Garnett_sp_eggs3.csv')[c("HANZAB.Name","taxon" )]
a<-merge(a,HANZAB,by.x='Var1',by.y='taxon',all.x=TRUE )
#write.csv(a,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/speciesOfInterest_',
#                   as.Date(Sys.time()),'.csv'),row.names=FALSE)

#write out observations
write.csv(eggdat, paste0('/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/UniqueLatLongYearMonthSourceEggobs_',
                         as.Date(Sys.time()),'.csv'),row.names=FALSE)







