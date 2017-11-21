rm(list = ls())

#prep species list
 
sp<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/speciesOfInterest_2015-07-06.csv')
colnames(sp)[1]<-"ScientificName"
colnames(sp)[4]<-"CommonName"
sp$genus<-t(as.data.frame(strsplit(as.character(sp$ScientificName),
                            " ")))[,1]
sp<-sp[,c("ScientificName","CommonName","genus","Family.scientific.name","Order","Freq","SpNo")]
colnames(sp)<-c("ScientificName","CommonName","genus","family","order","Freq","SpNo")

#page numbers
pg<-read.csv('/Users/daisy/Google Drive/PhD/birdTraits/IncubationFledging/speciesOfInterest.csv')[,c('Species','book','pgNu')]
sp<-merge(sp,pg,by.x="ScientificName",by.y="Species",all.x=TRUE)

#Garnett traits - 
gnt<-read.csv("/Users/daisy/Google Drive/PhD/birdTraits/Garnett paper/NSD-Data Descriptor/Australian_Bird_Data_Version_1_0.csv")
gnt$Taxon.scientific.name<- str_trim(gnt$Taxon.scientific.name)
gnt$Taxon.scientific.name<-capitalize(tolower(gnt$Taxon.scientific.name))
gnt<-subset(gnt,Species ==1)
gnt<- subset(gnt,Population.description == "Endemic (breeding only)" 
             | Population.description == "Australian" 
             | Population.description == "Endemic (entirely Australian)" 
             | Population.description =="Introduced" )

gnt<-gnt[,c("Taxon.scientific.name","Population.description","EPBC.Status.Nov..2014" ,
            "Australian.IUCN.Red.List.status.2014","Australian.IUCN.Red.List.criteria.2014",
            "Global.IUCN.status.2014","Global.IUCN.criteria.2014","Breeding.habitat..Arid.shrubland",                                           
            "Breeding.habitat..Chenopod.shrubland","Breeding.habitat..Heath","Breeding.habitat..Triodia.hummock.grassland",                                
            "Breeding.habitat..Other.grassland","Breeding.habitat..Mallee","Breeding.habitat..Tropical.savanna.woodland",                                
            "Breeding.habitat..Temperate.dry.scleorphyll.forest.and.woodland","Breeding.habitat..Temperate.wet.scleorphyll.forest.and.woodland",            
            "Breeding.habitat..Rainforest","Breeding.habitat..Mangrove","Breeding.habitat..inland.wetland",                                           
            "Breeding.habitat..Beaches.and.sand.cays","Breeding.habitat..Rocky.coasts.and.islets",                                 
            "Breeding.habitat..Other.non.Australian.habitats","Patch.clade", "Hackett.fine.clades","Hackett.coarse.clades")]

sp<-merge(sp,gnt,by.x="ScientificName",by.y="Taxon.scientific.name",all.x=TRUE)
#get traits of interest
traits<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/speciesOfInterest_2015-05-05.csv")[c("Species",
            "ClutchSizeMean","clutchDerived","RateOfLay","RateOfLayDerived","IncubationMean","IncubationDerivedHBW",
            "IncubationDerivedCloselyRelatedSpecies","FledgingMean","FledgingDerivedHBW",
            "FledgingDerivedCloselyRelatedSpecies")]

sp<-merge(sp,traits,by.x="ScientificName",by.y="Species",all.x=TRUE)


write.csv(sp,'/Users/daisy/Google Drive/PhD/BreedingTiming/tables/speciesOfInterest_2015-07-15.csv',row.names=FALSE)
