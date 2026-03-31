rm(list=ls())

#This code converts stroma to pixels to (x,y) coordinates
#Step 1: Plot stroma in BW, crop image, save image 
#Step 2: Convert BW image to pixels
#Step 3: Convert pixel back to coordinates for spatial statistics 
#Output with coordinates: "number_Converted_coordinates_DIAMETER_TRUE.csv"

#load libraries
library(sf)
library(dplyr)
library(ggplot2)
library(png)
library(magick)
#library(rstudioapi)
library(readr)

#current_path = rstudioapi::getActiveDocumentContext()$path 
#setwd(dirname(current_path ))
setwd(system("pwd", intern = T) )

myFile <- paste(readLines("fileName.txt"), collapse=" ")    
mynumber <- paste(readLines("number.txt"), collapse=" ") 
number <-as.numeric(mynumber)

#setup and import data
diameter <- TRUE #Whether to pixelate based on cell diameter 
#number <- 413910 #402016 #Describes number associated with experiment
scalefactor <- 0.5022 #micron per pixel
results <- readRDS(myFile)

name <- names(results)
regions <- results[[name]]$regions
ellipses<- results[[name]]$ellipses
regions_stroma<-filter(regions, class_label == "Stroma_fibr" )

ll = 1
for(yyq in regions_stroma$geometry){
  if ( is.na(sf::st_is_valid(yyq)) ){print(ll)
    sapply(regions_stroma$geometry[ll], function(x) nrow(x[[1]]))
    g = regions_stroma
    g$geometry[ll] = NULL
    regions_stroma = g}
  ll <- ll + 1
}
# sf::st_is_valid(regions_stroma) 
#sf::st_make_valid(regions_stroma$geometry[2146])
#sapply(regions_stroma$geometry[582], function(x) nrow(x[[1]]))
#g = regions_stroma
#g$geometry[582] = NULL
#regions_stroma = g
#plot(regions_stroma)

#parameters
cell_diam <- 15 #microns

#Create directories to save plots
# if (!dir.exists(paste(getwd(),"/results/",number,sep=""))){
#   dir.create(paste(getwd(),"/results/",number,sep=""))
# }
# if (!dir.exists(paste(getwd(),"/results/",number,"/whole_tissue/",sep=""))){
#   dir.create(paste(getwd(),"/results/",number,"/whole_tissue/",sep=""))
# }

#####################################################
# Step 1: Plot stroma in BW, crop image, save image #
#####################################################

#Plot stroma in BW
stroma_plot <- ggplot() +
  geom_sf(data = regions_stroma,
          fill="black",
          lwd=0)+
  theme_bw() + 
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.border=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),legend.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(x=NULL, y=NULL)
print(stroma_plot)

#Save plot
stroma_plot_filename <- paste(number,"_stroma_plot.png",sep="")
ggsave(stroma_plot_filename, stroma_plot) #,width=0.5,height=0.5,)

#Read plot
imread <- image_read(stroma_plot_filename)
imread <- image_trim(imread)

#Save cropped plot
cropped_stroma_plot_filename <- paste(number,"_cropped_stroma_plot.png",sep="")
image_write(imread,cropped_stroma_plot_filename,format="png")

#Read cropped plot
cropped_imread <- image_read(cropped_stroma_plot_filename) #produces raster image

######################################
# Step 2: Convert BW image to pixels #
######################################

#Determine bounding box coordinates of histology and convert to microns (scalefactor)
hist_minx <- as.numeric(st_bbox(regions_stroma)$xmin)*scalefactor
hist_maxx <- as.numeric(st_bbox(regions_stroma)$xmax)*scalefactor
hist_miny <- as.numeric(st_bbox(regions_stroma)$ymin)*scalefactor
hist_maxy <- as.numeric(st_bbox(regions_stroma)$ymax)*scalefactor

#resize image so that 1 px length = cell_diam microns
px_x <- round((hist_maxx-hist_minx)/cell_diam) #number of cells in x
px_y <- round((hist_maxy-hist_miny)/cell_diam) #number of cells in y
resized_matrix <- image_scale(cropped_imread,paste("",px_x,"",sep="")) #resize proportionally to width

#save small pixelated image, read image and convert to matrix of 0 and 1
pixel_stroma_plot_filename <- paste(number,"_pixel_stroma_plot.png",sep="")
if(diameter){
  image_write(resized_matrix,pixel_stroma_plot_filename,format="png")
  img <- readPNG(pixel_stroma_plot_filename) #produces raster image
  BWmatrix <- round(img) #0 if x<=0.5, 1 if x>0.5 #Threshold to define stroma (0) vs. empty (1)
} else {
  img <- readPNG(cropped_stroma_plot_filename)
  BWmatrix <- round(img) #0 if x<=0.5, 1 if x>0.5 #Threshold to define stroma (0) vs. empty (1)
}

Stroma <- which(apply(BWmatrix, c(1,2), function(x) x<=0.5),arr.ind=TRUE) 
Stroma <- Stroma[,ncol(Stroma):1] 
colnames(Stroma) <- c("X","Y")

###################################################################
#Step 3: Convert pixel back to coordinates for spatial statistics #
###################################################################

#Transform data so that coordinate is at center of pixel
Stroma[,"X"] <- sapply(Stroma[,"X"],function(x) hist_minx+(hist_maxx-hist_minx)/(2*px_x)+(x-1)*(hist_maxx-hist_minx)/(px_x))
Stroma[,"Y"] <- sapply(Stroma[,"Y"],function(x) hist_maxy-(hist_maxy-hist_miny)/(2*px_y)-(x-1)*(hist_maxy-hist_miny)/(px_y))

#Attach to cell coordinates in data frame
coordinates<-as.data.frame(Stroma)

#Attach other cell coordinates for plot
coordinates<-cbind(class_label="Stroma",coordinates)
cells<- as.data.frame(st_coordinates(ellipses)*scalefactor)
cells<-cbind(class_label=ellipses$class_label,cells)
coordinates<-rbind(coordinates,cells)

#csv file of cell coordinates (pixel->(x,y))
write.table(coordinates, file = paste(number,"_Converted_coordinates_DIAMETER_",diameter,".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

#Un-comment to show plot
#plot pixels on top of stroma
# p <- ggplot() +
#   geom_sf(data = regions_stroma,
#           aes(fill = class_label),
#           fill="blue",
#           alpha = 0.3,
#           #lwd = 0.5) +
#           lwd=0) +
#   geom_tile(data=subset(coordinates,class_label=="Stroma"),aes(x=X/scalefactor,y=Y/scalefactor),fill="blue")+#,alpha = 1)+# size=1)+ #+#shape="."
#   geom_sf(data = ellipses,
#           aes(color=class_label),
#           alpha = 0.5,size=0) + #, size = 0.1) +
#   geom_vline(xintercept = hist_minx/scalefactor)+
#   geom_vline(xintercept = hist_maxx/scalefactor)+
#   geom_hline(yintercept = hist_miny/scalefactor)+
#   geom_hline(yintercept = hist_maxy/scalefactor)+
#   theme(legend.position = 'bottom')+ #+
#   xlim(6000,8000)+
#   ylim(8000,10000)
# print(p)
# 
# ggsave(paste(getwd(),"/results/",number,"/whole_tissue/",number,"_allCells_DIAMETER_",diameter,".png",sep=""), p) #,width=0.5,height=0.5,)


