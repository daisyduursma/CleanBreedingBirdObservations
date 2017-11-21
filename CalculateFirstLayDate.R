

rm(list = ls())

library(Hmisc)
library(car)
library(zoo) 
library(stringr)
library(circular)

#return the julian date of observation
CalcEpoch <- function(PointLay){
  epoch<-as.numeric(strftime(as.Date(PointLay, origin = "1970-01-01"),format = "%j"))
  return(epoch)
}

#return the year of observation
CalcEpochYear <- function(PointLay){
  yr<-as.numeric(strftime(as.Date(PointLay, origin = "1970-01-01"),format = "%Y"))
  return(yr)
}
#return the month of observation
CalcEpochMonth <- function(PointLay){
  mn<-as.numeric(strftime(as.Date(PointLay, origin = "1970-01-01"),format = "%m"))
  return(mn)
}

CalcbinSize <- function(epochDates){
  d <- density(epochDates,na.rm = TRUE,from = 1, to = 366)
  ap <- approxfun(x=d$x, y=d$y)
  binSize<-integrate(ap, 1, 30.5)[[1]]
  for (mnth in 1:11){
    binSize<-c(binSize,integrate(ap, 30.5*mnth, (30.5*mnth)+30.5)[[1]])
  }
  return(binSize)
}

quantile_df <- function(x, probs, na.rm =F, names = F, type = 7, ...){
  z <- quantile(x, probs, na.rm, names, type)
  return(t(data.frame(row.names = probs, values = z)))
}

chi_sq <- function (x,y){
  #check vectorssame length and return e
  if(length(x) != length(y)) stop("x and y have different lengths")
  chi_val<-list()
  for (i in 1:length(x)) {
    chi_val[[i]]<-(y[i]-x[i])^2/x[i]
  }
   #remove Na and Inf
  chi<-na.omit(data.frame(do.call("rbind",chi_val))[1])[,1]
  chi[is.infinite(chi)] <- 0
  return(sum(chi))
}


dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
#TODO add new file
obs<-read.csv(paste0(dat.dir,'UniqueLatLongYearMonthSourceEggobs_2015-10-07.csv'))
obs$Scientific.Name<-as.character(obs$Scientific.Name)


#get list of incubation, fledging, clutch size etc

inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
obs<-obs[obs$Scientific.Name %nin% RMspecies,] #remove unwanted species


#species names
species<-as.character(inc$ScientificName)
#subset observations by species interested in
obs<-obs[obs$Scientific.Name %in% species,]



fld<-list()
for(i in 1:length(species)){
  spDat<-subset(obs,Scientific.Name==species[i])
  incub<-subset(inc,ScientificName==species[i])[1,'IncubationMean']
  fledge<-subset(inc,ScientificName==species[i])[1,'FledgingMean']
  #assume incubation starts on day last egg laid
  lay<-subset(inc,ScientificName==species[i])[1,'RateOfLay']*
          (subset(inc,ScientificName==species[i],select = 'ClutchSizeMean',drop=TRUE)-1)
  #first eggdate
  spDat$PLegg<-floor(apply(cbind(spDat$startEgg-lay,spDat$endEgg-incub-lay),
                           1, mean, trim = 0))
  spDat$PLstartEgg<-floor(spDat$startEgg-lay)
  spDat$PLendEgg<-floor(spDat$endEgg-incub-lay)
  #museuam egg
  spDat$PLeggMuseum_full<-floor(apply(cbind(spDat$startEgg-lay,spDat$endEgg-(incub)-lay), 1, mean, trim = 0))
  spDat$PLeggMuseum_fourth<-floor(apply(cbind(spDat$startEgg-lay,spDat$endEgg-(incub/4)-lay), 1, mean, trim = 0))
  #Hatch
  spDat$PLhatch<-spDat$startHatch-incub-lay
 #first young
  spDat$PLstartYoung<-spDat$startYoung-incub-lay
  spDat$PLstartYoung_fullFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-fledge), 1, mean, trim = 0))
  spDat$PLstartYoung_halfFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-(fledge/2)), 1, mean, trim = 0))
  spDat$PLstartYoung_fourthFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-(fledge/4)), 1, mean, trim = 0))
  #unknown
  spDat$PLUnkn_fullFlg<-floor(apply(cbind(spDat$startUnknown-lay,spDat$startUnknown-incub-lay,spDat$startUnknown-fledge-incub-lay), 1, mean, trim = 0))
  spDat$PLUnkn_halfFlg<-floor(apply(cbind(spDat$startUnknown-lay,spDat$startUnknown-incub-lay,spDat$startUnknown-(fledge/2)-incub-lay), 1, mean, trim = 0))
  spDat$PLUnkn_fourthFlg<-floor(apply(cbind(spDat$startUnknown-lay,spDat$startUnknown-incub-lay,spDat$startUnknown-(fledge/4)-incub-lay), 1, mean, trim = 0))
  spDat$PLUnkn_zeroFlg<-floor(apply(cbind(spDat$startUnknown-lay,spDat$startUnknown-incub-lay), 1, mean, trim = 0))
  fld[[i]]<-spDat
  message(i)
}
finObs<-do.call("rbind",fld)


###### make final point of lay for NRS data
###### uses information about nest build, timing of young, hatch
finObs$PL<-NA #make PL equal to NA
#find the point of lay between start of fledge, eggs and hatch, this is the mean between minimum difference of any of the backcalculated points
finObs$PL<-with(finObs,ifelse(sourceName=='NestRecordScheme',
                apply( matrix(c(PLhatch,PLstartEgg,PLendEgg,PLstartYoung,endBuild),ncol=5),1, #set up matrix
                      function(x) ( (sort(na.omit(x))[(which(abs(diff(sort(na.omit(x))))==min(abs(diff(sort(na.omit(x)))))))[1]] + sort(na.omit(x))[(which(abs(diff(sort(na.omit(x))))==min(abs(diff(sort(na.omit(x)))))))[1]+1]) /2) )
                  ,PL)) #if the above does not have values, NA
finObs$Acc<-with(finObs,ifelse(sourceName=='NestRecordScheme',
                          apply(matrix(c(PLhatch,PLstartEgg,PLendEgg,PLstartYoung,endBuild),ncol=5),1, 
                                function(x) min(abs(diff(sort(na.omit(x)))))/2),NA))
#fix obs where only hatch is given,PL is same as PLhatch and Acc is 0
finObs$Acc<-with(finObs,ifelse((sourceName=='NestRecordScheme'& is.na(PL) & !is.na(PLhatch)),0,Acc))
finObs$PL<-with(finObs,ifelse((sourceName=='NestRecordScheme'& is.na(PL) & !is.na(PLhatch)),PLhatch,PL))

                              
#make inf set to NA and round 
finObs$Acc[is.infinite(finObs$Acc)] <- NA
finObs$Acc<-round(finObs$Acc)
#if startegg and endegg are within 2 day, Acc > 10 days and no other info is given, convert to solitary egg.
finObs$eggtime<-finObs$startEgg-finObs$endEgg
finObs$type<- with(finObs,ifelse((!is.na(Acc) & Acc > 10 & eggtime < 3 & sourceName=='NestRecordScheme'),"solitary-egg",type))
finObs$Acc<- with(finObs,ifelse((!is.na(Acc) & Acc > 10 & eggtime < 3 & sourceName=='NestRecordScheme'),NA,Acc))

#make sure ones with hatch only have accuracy of 0 and are coded as multi-obs
finObs$type<- with(finObs,ifelse((!is.na(Acc) & Acc <= 10 & sourceName=='NestRecordScheme' &  type =='multi-egg'),"multi-obs",type))
finObs$type<- with(finObs,ifelse((!is.na(Acc) & Acc <= 10 & sourceName=='NestRecordScheme' &  type == "youngOut"),"multi-obs",type))
#make sure single egg with build are multi-obs
finObs$type<- with(finObs,ifelse((!is.na(Acc) & Acc <= 10 & sourceName=='NestRecordScheme' &  !is.na(endBuild) & !is.na (startEgg))
                                 ,"multi-obs",type))
#check to see that last build date is not later than Point of Lay
finObs$PL<-with(finObs,ifelse((!is.na(endBuild) & !is.na(PL) & (endBuild > PL) & sourceName=='NestRecordScheme'),
      endBuild,PL))



######Convert the year, month, day since 1970 to day of the year
finObs$DOY_PL<-CalcEpoch(finObs$PL)
finObs$DOY_PLeggMuseum_full<- CalcEpoch(finObs$PLeggMuseum_full)
finObs$DOY_PLeggMuseum_fourth<-CalcEpoch(finObs$PLeggMuseum_fourth)
finObs$DOY_youngfull<-CalcEpoch(finObs$PLstartYoung_fullFlg)
finObs$DOY_younghalf<-CalcEpoch(finObs$PLstartYoung_halfFlg)
finObs$DOY_youngfourth<-CalcEpoch(finObs$PLstartYoung_fourthFlg)
finObs$DOY_youngzero<-CalcEpoch(finObs$PLstartYoung)
finObs$DOY_PLUnknFull<-CalcEpoch(finObs$PLUnkn_fullFlg)
finObs$DOY_PLUnknhalf<-CalcEpoch(finObs$PLUnkn_halfFlg)
finObs$DOY_PLUnknfourth<-CalcEpoch(finObs$PLUnkn_fourthFlg)
finObs$DOY_PLUnknzero<-CalcEpoch(finObs$PLUnkn_zeroFlg)

#####subset data for fitting density and calculating chi-square
finObs$ID<-1:nrow(finObs) #make ID for subsettting
NRSObsAcc10<-subset(finObs, sourceName=='NestRecordScheme' & type=='multi-obs' & Acc<=10)
NRSAccObs<-subset(finObs, sourceName=='NestRecordScheme' & Acc<=5 & type=='multi-obs')
finObsOther<-finObs[finObs$ID %nin% NRSObsAcc10$ID,] #all other observations except NRS multi

#########################
calcChiSquare<-TRUE

if (calcChiSquare == TRUE) {
  #empty list
  alldat <- list()
  #for species
  for (sp in 1:length(species)) {
    #Get species that are of interest
    #for Species, Genus, Family, Order, or  Altricial find best method to convert to Point of Lay
    AP <- as.vector(subset(inc,ScientificName == species[sp])$DoD)
    Ord <- as.vector(subset(inc,ScientificName == species[sp])$order)
    Fam <- as.vector(subset(inc,ScientificName == species[sp])$family)
    gen <- as.vector(subset(inc,ScientificName == species[sp])$genus)
    
    #test for species species
    Dodobs <-
      finObsOther[finObsOther$Scientific.Name %in% species[sp],]
    #Get NRS accurate, observations for those species
    NRSsub <- NRSAccObs[NRSAccObs$Scientific.Name %in% species[sp],]
    NRS <- NRSsub$DOY_PL
    NRScount <- nrow(NRSsub)
    ChiGroup <- paste('species',species[sp],sep = " - ")
    #test for genus
    if (NRScount < 100) {
      subspec <- as.character(subset(inc,genus == gen)[,'ScientificName'])
      Dodobs <-
        finObsOther[finObsOther$Scientific.Name %in% subspec,]
      NRSsub <- NRSAccObs[NRSAccObs$Scientific.Name %in% subspec,]
      NRS <- NRSsub$DOY_PL
      NRScount <- nrow(NRSsub)
      ChiGroup <- paste('genus',gen,sep = " - ")
    }
    if (NRScount < 100) {
      subspec <- as.character(subset(inc,family == Fam)[,'ScientificName'])
      Dodobs <-
        finObsOther[finObsOther$Scientific.Name %in% subspec,]
      NRSsub <- NRSAccObs[NRSAccObs$Scientific.Name %in% subspec,]
      NRS <- NRSsub$DOY_PL
      NRScount <- nrow(NRSsub)
      ChiGroup <- paste('family',Fam, sep = " - ")
    }
    if (NRScount < 100) {
      subspec <- as.character(subset(inc,order == Ord)[,'ScientificName'])
      Dodobs <-
        finObsOther[finObsOther$Scientific.Name %in% subspec,]
      NRSsub <- NRSAccObs[NRSAccObs$Scientific.Name %in% subspec,]
      NRS <- NRSsub$DOY_PL
      NRScount <- nrow(NRSsub)
      ChiGroup <- paste('order',Ord,sep = " - ")
    }
    if (NRScount < 100) {
      subspec <- as.character(subset(inc,DoD == AP)[,'ScientificName'])
      Dodobs <-
        finObsOther[finObsOther$Scientific.Name %in% subspec,]
      NRSsub <- NRSAccObs[NRSAccObs$Scientific.Name %in% subspec,]
      NRS <- NRSsub$DOY_PL
      NRScount <- nrow(NRSsub)
      ChiGroup <- paste('DegreeOfDev',AP,sep = " - ")
    }
    
    
    
    #get the back calcuated Points of lays
    DOY_PLeggMuseumFull <- na.omit(Dodobs[,'DOY_PLeggMuseum_full'])
    #DOY_PLeggMuseumhalf<-na.omit(Dodobs[,'DOY_PLeggMuseum_half'])
    #DOY_PLeggMuseumthird<-na.omit(Dodobs[,'DOY_PLeggMuseum_third'])
    DOY_PLeggMuseumfourth <-
      na.omit(Dodobs[,'DOY_PLeggMuseum_fourth'])
    DOY_youngfull <- na.omit(Dodobs[,'DOY_youngfull'])
    DOY_younghalf <- na.omit(Dodobs[,'DOY_younghalf'])
    #DOY_youngthird<-na.omit(Dodobs[,'DOY_youngthird'])
    DOY_youngfourth <- na.omit(Dodobs[,'DOY_youngfourth'])
    DOY_youngzero <- na.omit(Dodobs[,'DOY_youngzero'])
    DOY_PLUnknFull <- na.omit(Dodobs[,'DOY_PLUnknFull'])
    DOY_PLUnknhalf <- na.omit(Dodobs[,'DOY_PLUnknhalf'])
    #DOY_PLUnknthird<-na.omit(Dodobs[,'DOY_PLUnknthird'])
    DOY_PLUnknfourth <- na.omit(Dodobs[,'DOY_PLUnknfourth'])
    DOY_PLUnknzero <- na.omit(Dodobs[,'DOY_PLUnknzero'])
    
          ############density curves of point of lay
           {
            pdf(paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/figures/FittingPointOfLay/FittingPointOfLay_',species[sp],'_',as.Date(Sys.Date()),'.pdf'), 8,5)
            par(mar=c(0,2,3,1))
            layout(matrix(c(1,1,2,3,4,5,6,6),4,2,byrow = TRUE),heights = c(1,4,4,2))
            #upper label
            plot(1, type = "n", axes=FALSE, xlab="", ylab="", main = paste0(species[sp],", ",ChiGroup))
            #NRS
            par(mar=c(4,4,2,2))
            d <- density(NRS)
            plot(d, type="n", main="Nest Record Scheme",xlab="day",xaxt='n',xlim=c(1,366))
            polygon(d, col="grey", border="grey")
            #hist(NRS,breaks=seq(0,366,by=6), main ="Nest Record Scheme",xlab="day",xaxt='n',xlim=c(0,365),prob=TRUE)
            #lines(density(NRS), col="blue", lwd=2,lty=2)
            #Museum Egg
            plot(d, type="n", main ="Solitary eggs",xlab="day",xaxt='n',xlim=c(1,366))
            polygon(d, col="grey", border="grey")
            #hist(NRS,breaks=seq(0,366,by=6), main ="Museum eggs",xlab="day",xaxt='n',xlim=c(1,365),prob=TRUE)
            #lines(density(NRS), col="blue", lwd=2,lty=2)
            lines(density(DOY_PLeggMuseumFull), col="green", lwd=2)
            #lines(density(DOY_PLeggMuseumhalf), col="deepskyblue", lwd=2)
            #lines(density(DOY_PLeggMuseumthird), col="deepskyblue4", lwd=2)
            lines(density(DOY_PLeggMuseumfourth), col="purple3", lwd=2)
            #NRS Young
            plot(d, type="n", main ="Young",xlab="day",xaxt='n',xlim=c(1,366))
            polygon(d, col="grey", border="grey")
            #   hist(NRS,breaks=seq(0,366,by=6), main ="NRS young",xlab="day",xaxt='n',xlim=c(0,365),prob=TRUE)
            #   lines(density(NRS), col="blue", lwd=2,lty=2)
            lines(density(DOY_youngfull), col="green", lwd=2)
            lines(density(DOY_younghalf), col="deepskyblue", lwd=2)
            #lines(density(DOY_youngthird), col="deepskyblue4", lwd=2)
            lines(density(DOY_youngfourth), col="purple3", lwd=2)
            lines(density(DOY_youngzero), col="black", lwd=2)
            #ATLAS
            plot(d, type="n", main ="Unknown",xlab="day",xaxt='n',xlim=c(1,366))
            polygon(d, col="grey", border="grey")
            #   hist(NRS,breaks=seq(0,366,by=6), main ="ATLAS breeding",xlab="day",xaxt='n',xlim=c(0,365),prob=TRUE)
            #   lines(density(NRS), col="blue", lwd=2,lty=2)
            lines(density(DOY_PLUnknFull), col='green', lwd=2)
            lines(density(DOY_PLUnknhalf), col='deepskyblue', lwd=2)
            #lines(density(DOY_PLUnknthird), col='deepskyblue4', lwd=2) #these have a better fit
            lines(density(DOY_PLUnknfourth), col='purple3', lwd=2) #need to decide before third and fourth
            lines(density(DOY_PLUnknzero), col='black', lwd=2)
            #legend
            par(mar=c(0,2,3,1))
            plot(1, type = "n", axes=FALSE, xlab="", ylab="")
            legend(x = "top",inset = 0,
                   legend = c("NRS","full","half","fourth","zero"),
                   col=c("grey","green", "deepskyblue","purple3","black"),
                   lwd=5, cex=.9, horiz = TRUE)
            dev.off()
          }
    
    qdf <-
      as.data.frame(t(c(
        species[sp], ChiGroup, NRScount,nrow(Dodobs)
      )))
    colnames(qdf) <-
      c("Species","CorrectionLevel","NRSobs","OtherObs")
    #     #calculate size of 12 (monthly) bins
    binsDF <- list()#make list to feed them to
    binsDF[[1]] <- CalcbinSize(NRS)
    binsDF[[2]] <- CalcbinSize(DOY_PLeggMuseumFull)
    binsDF[[3]] <- CalcbinSize(DOY_PLeggMuseumfourth)
    binsDF[[4]] <- CalcbinSize(DOY_youngfull)
    binsDF[[5]] <- CalcbinSize(DOY_younghalf)
    binsDF[[6]] <- CalcbinSize(DOY_youngfourth)
    binsDF[[7]] <- CalcbinSize(DOY_youngzero)
    binsDF[[8]] <- CalcbinSize(DOY_PLUnknFull)
    binsDF[[9]] <- CalcbinSize(DOY_PLUnknhalf)
    binsDF[[10]] <- CalcbinSize(DOY_PLUnknfourth)
    binsDF[[11]] <- CalcbinSize(DOY_PLUnknzero)
    MonthDat <-
      as.data.frame(do.call("rbind",binsDF),row.names = '')#put them together as dataframe
    colnames(MonthDat) <- paste0("Month",1:12)#add column names
    qdf <- cbind(qdf,MonthDat)#add to summary table
    qdf$scenario <-
      c(
        'NRS','MusEggFullIncubation','MusEggfourthIncIncubation','NRS_youngFull_fledge',
        'NRS_younghalf_fledge','NRS_youngfourth_fledge','NRS_youngzero_fledge',
        'ATLAS_UnknFull_fledge','ATLAS_Unknhalf_fledge','ATLAS_Unknfourth_fledge','ATLAS_Unknzero_fledge'
      )
    
    #http://www.cliffsnotes.com/math/statistics/bivariate-relationships/chi-square-x2
    #calc Chi-square value between expected and observed.
    #This is not a test but rather a ranking
    qdf$chi_sq <- NA
    for (qdfRow in 1:nrow(qdf))
    {
      x <-
        as.numeric(subset(qdf,scenario == 'NRS',select = paste0("Month",1:12)))
      y <-
        as.numeric(subset(qdf,select = paste0("Month",1:12))[qdfRow,])
      qdf$chi_sq[qdfRow] <- chi_sq(x,y)
    }
    
    
    alldat[[sp]] <- qdf
    
    message(sp)
    
  }
  finDat <- as.data.frame(do.call("rbind",alldat))

  #get amounts to use
  proportions<-list()
  for(sp in 1:length(species)){
    subfam<-subset(finDat,Species==species[sp])
    CorrectionLevel<-as.vector(unique(subfam$CorrectionLevel))
    egg<-subset(subfam, chi_sq==min(subset(subfam, scenario == "MusEggFullIncubation" |
               scenario == "MusEggfourthIncIncubation")$chi_sq), scenario)[1,]
    
    young<-subset(subfam, chi_sq==min(subset(subfam, scenario == "NRS_youngFull_fledge" |
                                                scenario == "NRS_younghalf_fledge"|
                                                  scenario == "NRS_youngfourth_fledge"|
                                                  scenario == "NRS_youngzero_fledge")$chi_sq), scenario)[1,]
    
    unknown<-subset(subfam, chi_sq==min(subset(subfam, scenario == "ATLAS_UnknFull_fledge" |
                                                   scenario == "ATLAS_Unknhalf_fledge"|
                                                   scenario == "ATLAS_Unknfourth_fledge"|
                                                   scenario == "ATLAS_Unknzero_fledge")$chi_sq), scenario)[1,]
  
    
    proportions[[sp]]<-cbind(species[sp],CorrectionLevel,egg,young,unknown)
  
    message(sp)
    }
  finprop<-as.data.frame(do.call("rbind",proportions))
  colnames(finprop)[1]<-'Species'
  
  
  #write.csv(finDat,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/QuantilesChiSquareBreedingDataTypes_family_10dayAccuracy',as.Date(Sys.Date()),'.csv'),row.names=FALSE)
  write.csv(finprop,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/Proportions_5day',as.Date(Sys.Date()),'.csv'),row.names=FALSE)
}


##########Based on Chi-square values make final DOY_PL values 
#finObsOther - needs DOY_PL and NRSMulti does not

#PLvals<-read.csv("//Users/daisy/Google Drive/PhD/BreedingTiming/tables/Proportions2015-05-26.csv")
PLvals<-finprop

#get list of Orders
#family<-unique(inc$Family.scientific.name.y) #get list of Orders

finPL<-list()
for (sp in 1:length(species)){
  #spec<-species[sp]
  #Get observations for those species
  OrdObs<-finObsOther[finObsOther$Scientific.Name %in% species[sp],] 
  vals<-subset(PLvals,Species == paste(species[sp]))
  #eggs
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$egg =="MusEggfourthIncIncubation" & !is.na(DOY_PLeggMuseum_fourth)),DOY_PLeggMuseum_fourth,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$egg =="MusEggFullIncubation" & !is.na(DOY_PLeggMuseum_full)),DOY_PLeggMuseum_full,DOY_PL))#calculated difference between hatch and egg dated
  #young
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$young =="NRS_youngzero_fledge" & !is.na(DOY_youngzero)),DOY_youngzero,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$young =="NRS_youngfourth_fledge" & !is.na(DOY_youngfourth)),DOY_youngfourth,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$young =="NRS_younghalf_fledge" & !is.na(DOY_younghalf)),DOY_younghalf,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$young =="NRS_youngFull_fledge" & !is.na(DOY_youngfull)),DOY_youngfull,DOY_PL))#calculated difference between hatch and egg dated
  #unknown
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$unknown =="ATLAS_Unknzero_fledge" & !is.na(DOY_PLUnknzero)),DOY_PLUnknzero,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$unknown =="ATLAS_Unknfourth_fledge" & !is.na(DOY_PLUnknfourth)),DOY_PLUnknfourth,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$unknown =="ATLAS_Unknhalf_fledge" & !is.na(DOY_PLUnknhalf)),DOY_PLUnknhalf,DOY_PL))#calculated difference between hatch and egg dated
  OrdObs$DOY_PL<-with(OrdObs,ifelse((vals$unknown =="ATLAS_UnknFull_fledge" & !is.na(DOY_PLUnknFull)),DOY_PLUnknFull,DOY_PL))#calculated difference between hatch and egg dated

  finPL[[sp]]<-OrdObs
}

finPL<-as.data.frame(do.call("rbind",finPL))

#get years
finPL$year<-CalcEpochYear(finPL$PLegg)
finPL$month<-CalcEpochMonth(finPL$PLegg)
finPL$year<-with(finPL,ifelse(is.na(year),CalcEpochYear(finPL$startYoung),year))
finPL$month<-with(finPL,ifelse(is.na(month),CalcEpochMonth(finPL$startYoung),month))
finPL$year<-with(finPL,ifelse(is.na(year),CalcEpochYear(finPL$startUnknown),year))
finPL$month<-with(finPL,ifelse(is.na(month),CalcEpochMonth(finPL$startUnknown),month))
NRSObsAcc10$year<-CalcEpochYear(NRSObsAcc10$PL)
NRSObsAcc10$month<-CalcEpochMonth(NRSObsAcc10$PL)

finPL<-rbind(finPL,NRSObsAcc10)


#remove DOY_PL where there is only NRS Build data
table(is.na( finPL$DOY_PL))
finPL<-subset(finPL,!is.na(DOY_PL))

write.csv(finPL,paste0(dat.dir,'PointOfLayDayOfYear',
                       as.Date(Sys.time()),'.csv'),row.names=FALSE)








##################################Results in Paper

#working directory
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
finPL<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2015-07-09.csv'))

#make figure of density of observations in years for 3 types of data
NRSYear <-subset(finPL, sourceName == "NestRecordScheme")$year
musYear<-subset(finPL, sourceName == "SouthAustraliaMuseaum" |
              sourceName == "QueenVictoria"| 
              sourceName ==  "AustraliaMuseum"|
              sourceName == "MuseumVictoria" |
              sourceName == "AustralianNationalWildlifeCollection" |
              sourceName == "WesternAustraliaMuseum"|
              sourceName == "TasmanianMuseumArtGallery" |
              sourceName == "NorthernTerritoryMuseumAndArtGallery" |
              sourceName == "QueenslandMuseum") $year 
altYear<-finPL[grep("^ATLAS", finPL$sourceName), ]$year
ABBBSyear<-subset(finPL, sourceName == "ABBBS")$year
ebirdyear<-subset(finPL, sourceName == "eBIRD")$year
plot(density(altYear), main ="PDF year of collection",xlab="year",xlim=c(1770,2015),col="purple3", lwd=2)
lines(density(NRSYear), col="green", lwd=2)
lines(density(musYear), col="deepskyblue", lwd=2)
lines(density(ABBBSyear), col="darkred", lwd=2)
lines(density(ebirdyear), col="darkgoldenrod2", lwd=2)
legend("topleft",c("Atlas","NRS","Museum","ABBBS","eBIRD"),text.col=c("purple3","green","deepskyblue","darkred","darkgoldenrod2"))

###################################
#get mean species observation rate by order
aggregate(Freq ~ Order.y, data = inc2,mean)
#get mean species observation rate by family
aggregate(Freq ~ Family.scientific.name.y, data = inc2,mean)
#mean and sd across all species
mean(inc2$Freq)
sd(inc2$Freq)

#find orders that do not have at least 100 observations in NRS accurate data
DoDev<-unique(inc$Order.y)
NRSCount<-list()
for (dev in 1:length(DoDev)){
  DoDspec<-as.character(subset(inc,Order.y == DoDev[dev])[,'Species'])
  NRSsub<-NRSAccObs[NRSAccObs$Scientific.Name %in% DoDspec,] 
  NRS<-NRSsub$DOY_PL
  NRSCount[[dev]]<-nrow(NRSsub)
}  
NRSCount<-cbind(as.data.frame(DoDev),do.call("rbind",NRSCount))

  

#differnces in the number of days from using differnt proportions of incubation
mus<-subset(finPL, sourceName == "SouthAustraliaMuseaum" |
                  sourceName == "QueenVictoria"| 
                  sourceName ==  "AustraliaMuseum"|
                  sourceName == "MuseumVictoria" |
                  sourceName == "AustralianNationalWildlifeCollection" |
                  sourceName == "WesternAustraliaMuseum"|
                  sourceName == "TasmanianMuseumArtGallery" |
                  sourceName == "NorthernTerritoryMuseumAndArtGallery" |
                  sourceName == "QueenslandMuseum")

mean(mus$PLeggMuseum_full-mus$PLeggMuseum_fourth)
mean(mus$PLeggMuseum_fourth-mus$PLeggMuseum_half)
mean(mus$PLeggMuseum_fourth-mus$PLeggMuseum_third)

#differnces in the number of days from using differnt proportions of fledge
atl<-finPL[grep("^ATLAS", finPL$sourceName), ]
you<-atl<-finPL[grep("^ATLAS", finPL$sourceName), ]


mean(alt$PLUnkn_fullFlg-alt$PLUnkn_zeroFlg)
mean(alt$PLUnkn_fullFlg-alt$PLUnkn_halfFlg)
mean(alt$PLUnkn_fullFlg-alt$PLUnkn_thirdFlg)
mean(alt$PLUnkn_fullFlg-alt$PLUnkn_fourthFlg)





#first young
spDat$PLstartYoung<-spDat$startYoung-incub-lay
spDat$PLEggYoungDiff<-floor(apply(cbind(spDat$PLegg,spDat$PLstartYoung), 1, diff, trim = 0))
spDat$PLstartYoung_fullFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-fledge), 1, mean, trim = 0))
spDat$PLstartYoung_halfFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-(fledge/2)), 1, mean, trim = 0))
spDat$PLstartYoung_fourthFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-(fledge/4)), 1, mean, trim = 0))
spDat$PLstartYoung_thirdFlg<-floor(apply(cbind(spDat$PLstartYoung,spDat$PLstartYoung-(fledge/3)), 1, mean, trim = 0))





