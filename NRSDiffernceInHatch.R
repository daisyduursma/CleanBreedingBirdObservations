


#Look at variation in hatch date for differnt eggs


rm(list = ls())

library(raster)
library(car) 
library(stringr)
library(Hmisc)
library(lubridate)
library(gtools)

indir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/egg'


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
a<-as.data.frame(table(H$ID))

subset(a,Freq>1)

