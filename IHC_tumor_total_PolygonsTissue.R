rm(list=ls())
library(sp)
library(rmapshaper)
library(ggplot2)
library(magrittr)
library(sf)
library(spatstat)
library(goftest)
library(readr)
library(kSamples)
library(dplyr)
#library(rstudioapi)
#library(readr)


#current_path = rstudioapi::getActiveDocumentContext()$path 
#current_path <- getwd()
#setwd(dirname(current_path))
#setwd(getSrcDirectory()[1])
setwd(system("pwd", intern = T) )

myFile <- paste(readLines("fileName.txt"), collapse=" ")    
mynumber <- paste(readLines("number.txt"), collapse=" ") 
number <-as.numeric(mynumber)

diameter <- TRUE 
#scalefactor <- 0.5022 #micron per pixel
results <- readRDS(myFile)
name <- names(results)
regions <- results[[name]]$regions
#mydata <- read.csv(paste(getwd(),"/results/",number,"/whole_tissue/",number,"_Converted_coordinates_DIAMETER_",diameter,".csv",sep=""))

#Extract vertices of tumor_total 
#L1 refers to the main ring or holes, L2 to the ring id 
#(I think L1 defines unique polygons)
#tumor_total_owin <- as.data.frame(subset(regions,class_label=="IHC_tumor_total") %>% sf::st_coordinates())
tumor_total_owin <- as.data.frame(subset(regions,class_label=="IHC_tumor_total"))# %>% sf::st_coordinates())

polyList <- as.data.frame(tumor_total_owin$geometry)
myArea <- c()
for (t in polyList$geometry){myArea <- c(myArea, st_area(t))}
myMax <- max(myArea)
mv_simpl <- as.data.frame(st_simplify(subset(polyList$geometry, st_area(polyList$geometry) > myMax*0.05), 
                        preserveTopology = FALSE, dTolerance = 160)%>% sf::st_coordinates())
#plot(mv_simpl)
#Create list of lists of polygon vertices
all_list <- c()

# for(i1 in unique(tumor_total_owin[,"L1"])){
#   tumor_total_index <- subset(tumor_total_owin,L1==i1)
#   for (i2 in unique(tumor_total_index[,"L2"])){
#     poly_list <- list(list(x=subset(tumor_total_index,L2==i2)[nrow(subset(tumor_total_index,L2==i2)):1,"X"]*scalefactor, 
#                            y=subset(tumor_total_index,L2==i2)[nrow(subset(tumor_total_index,L2==i2)):1,"Y"]*scalefactor))
#     all_list <- rbind(all_list,poly_list)
#   }
# }
#write.csv(tumor_total_owin, file =( "Polygons.csv"))
write.csv(mv_simpl, file =( "PolygonsMainSmooth.csv"))

