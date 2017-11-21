#make supplementary tables for paper
dat.dir<-'/Users/daisy/Google Drive/PhD/BreedingTiming/tables'

breeding<-read.csv(paste0(dat.dir,'/BreedingQuantiles2015-06-02.csv'))
VanDerWal<-read.csv(paste0(dat.dir,'/VanDerWalSupT2.csv'))
VanDerWal$Scientific.Name<-paste(VanDerWal$genus,VanDerWal$species,sep=' ')
Garnett<-read.csv('/Users/daisy/Google Drive/PhD/birdTraits/Garnett paper/NSD-Data Descriptor/Australian_Bird_Data_Version_1_0.csv')

species<-breeding$Species 

#subset data for species
SubVan<<-VanDerWal[VanDerWal$Scientific.Name %in% species,c("Scientific.Name","Distance.km.","Direction.degrees.","Velocity.km.yr.")]
SubGar<<-Garnett[Garnett$Taxon.scientific.name %in% species,c("Taxon.scientific.name","Climate.specialisation","Specialisation.index","Sensitivity.index")]

findat<-merge(SubGar,SubVan,by.x="Taxon.scientific.name", by.y = "Scientific.Name",all.x=TRUE)

findat<-merge(breeding,findat,by.y="Taxon.scientific.name", by.x = "Species",all.x=TRUE)
findat$breedingDaysPer90<-NA

for (i in 1:nrow(findat)){
  startDate<-findat$Aus5[i]
  endDate<-findat$Aus95[i]
  
  findat$breedingDaysPer90[i]<-ifelse (endDate < startDate, 365-startDate+endDate,endDate - startDate )
  
}

write.csv(findat,paste0(dat.dir,"SummaryBreedingLengthSensitivity.csv"),row.names=FALSE)


plot(density(findat$breedingDaysPer90))


table(findat$breedingDaysPer90)
