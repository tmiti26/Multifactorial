rm(list=ls())
library(sp)
library(rmapshaper)
#library(ggplot2)
library(magrittr)
library(sf)
library(spatstat)
library(goftest)
library(readr)
library(kSamples)

#setup and import data
setwd(system("pwd", intern = T) )
myFile <- paste(readLines("fileName.txt"), collapse=" ")    
mynumber <- paste(readLines("number.txt"), collapse=" ") 
number <-as.numeric(mynumber)
results <- readRDS(myFile)


linn <-readLines(file("Rparams.txt",open="r"))
close(file("Rparams.txt",open="r"))
quadratSize <- as.numeric(linn[1]) #1000
PCF_r <- as.numeric(linn[2]) # 450
J_r <- as.numeric(linn[3]) #250
Z_r_length <- as.numeric(linn[4]) #500
dens_r_length <- as.numeric(linn[5]) #500
PCF_r_length <- as.integer(PCF_r*2)
J_r_length <- as.integer(J_r*2)

#Parameters
cell_diam <- 15 #13(Breast) #15(Lung) #microns
diameter <- TRUE #Whether to pixelate based on cell diameter 
scalefactor <- 0.5022 #micron per pixel

name <- names(results)
regions <- results[[name]]$regions
mydata <- read.csv(paste(getwd(),"/",number,"_Converted_coordinates_DIAMETER_",diameter,".csv",sep=""))

mysim <- function(ppp1,ppp2){
  superimpose(Pos = split(rlabel(ppp1))$Pos, Str = ppp2)
}

mysimN <- function(ppp1,ppp2){
  superimpose(Neg = split(rlabel(ppp1))$Neg, Str = ppp2)
}
## get the individual large polygons as analysis windows
tumor_total_owin_gem <- as.data.frame(subset(regions,class_label=="IHC_tumor_total"))
polyList <- as.data.frame(tumor_total_owin_gem$geometry)
myArea <- c()
for (t in polyList$geometry){myArea <- c(myArea, st_area(t))}
myMax <- max(myArea)
mv_simpl <- st_simplify(subset(polyList$geometry, st_area(polyList$geometry) > myMax*0.15), preserveTopology = FALSE, dTolerance = 180)
#plot(mv_simpl)
#No changes
averObs <- rep(0, PCF_r_length)
averMean <- rep(0, PCF_r_length)
averR <- rep(0, PCF_r_length)
averHi <- rep(0, PCF_r_length)
averLo <- rep(0, PCF_r_length)
averObsP <- rep(0, PCF_r_length)
averMeanP <- rep(0, PCF_r_length)
averRP <- rep(0, PCF_r_length)
averHiP <- rep(0, PCF_r_length)
averLoP <- rep(0, PCF_r_length)

rPCF <-seq(from = 0, to = PCF_r, by = 0.5)

averObsL <- rep(0, PCF_r_length)
averMeanL <- rep(0, PCF_r_length)
averRL <- rep(0, PCF_r_length)
averHiL <- rep(0, PCF_r_length)
averLoL <- rep(0, PCF_r_length)
averObsPL <- rep(0, PCF_r_length)
averMeanPL <- rep(0, PCF_r_length)
averRPL <- rep(0, PCF_r_length)
averHiPL <- rep(0, PCF_r_length)
averLoPL <- rep(0, PCF_r_length)

averObsJ <- rep(0, J_r_length)
averMeanJ <- rep(0, J_r_length)
averRJ <- rep(0, J_r_length)
averHiJ <- rep(0, J_r_length)
averLoJ <- rep(0, J_r_length)
averObsPJ <- rep(0, J_r_length)
averMeanPJ <- rep(0, J_r_length)
averRPJ <- rep(0, J_r_length)
averHiPJ <- rep(0, J_r_length)
averLoPJ <- rep(0, J_r_length)

rJ <-seq(from = 0, to = J_r, by = 0.5)

averObsG <- rep(0, J_r_length)
averMeanG <- rep(0, J_r_length)
averRG <- rep(0, J_r_length)
averHiG <- rep(0, J_r_length)
averLoG <- rep(0, J_r_length)
averObsPG <- rep(0, J_r_length)
averMeanPG <- rep(0, J_r_length)
averRPG <- rep(0, J_r_length)
averHiPG <- rep(0, J_r_length)
averLoPG <- rep(0, J_r_length)

rhoAver <- rep(0,Z_r_length)
rhoZ <- rep(0,Z_r_length)
rhoHi <- rep(0,Z_r_length)
rhoLo <- rep(0,Z_r_length)
rhoMean <- rep(0,Z_r_length)
rhoAverP <- rep(0,Z_r_length)
rhoZP <- rep(0,Z_r_length)
rhoHiP <- rep(0,Z_r_length)
rhoLoP <- rep(0,Z_r_length)
rhoMeanP <- rep(0,Z_r_length)

nncrossp <- c()
nncrossn <- c()

numb <- 0
for (i in 1:length(mv_simpl)){
  all_owin <- as.owin(mv_simpl[i]*scalefactor)
  #plot(all_owin)
  
  x <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('BrdU_pos')))]
  y <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('BrdU_pos')))]
  xy_inside<-inside.owin(x,y,all_owin)
  xx <-x[xy_inside]
  yy <- y[xy_inside]
  mypatternPos <- ppp(xx,yy,window=all_owin)#,marks=as.factor(mydata$class_label[xy_inside])) 
  
  xn <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('BrdU_neg')))]
  yn <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('BrdU_neg')))]
  xyn_inside<-inside.owin(xn,yn,all_owin)
  xxn <-xn[xyn_inside]
  yyn <- yn[xyn_inside]
  mypatternNeg <- ppp(xxn,yyn,window=all_owin)#,marks=as.factor(mydata$class_label[xy_inside])) 
  
  xs <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  ys <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  xys_inside<-inside.owin(xs,ys,all_owin)
  xxs <-xs[xys_inside]
  yys <- ys[xys_inside]
  mypatternStr <- ppp(xxs,yys,window=all_owin)#,marks=as.factor(mydata$class_label[xy_inside])) 
  
  print(mypatternPos$n/(mypatternNeg$n+mypatternPos$n))
  prolif <- mypatternPos$n/(mypatternNeg$n+mypatternPos$n)
  write.table(prolif, file = paste(i,"_myProlifFract",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
  
  
  #gettinf the edge of the stroma 
  vx <- c()
  vy <- c()
  indexStr <- unique(nncross(mypatternNeg, mypatternStr, k =2)$which)
  for (j in indexStr){
    vx <- c(vx, mypatternStr$x[j])
    vy <- c(vy, mypatternStr$y[j])}
  
  indexStrp <- unique(nncross(mypatternPos, mypatternStr, k = 2)$which)
  for (h in indexStrp){
    vx <- c(vx, mypatternStr$x[h])
    vy <- c(vy, mypatternStr$y[h])}
  
  mypatternStrEdge <- ppp(vx, vy,window= all_owin )
  
  # check to not analize too small off polygons
  if((mypatternStr$n + mypatternNeg$n + mypatternPos$n) > 700){
    numb <- numb + 1
    
    # #rho hat ##################################################
    Z <- distmap(mypatternStr)
    myRho <- rhohat(mypatternNeg,Z)
    myRhoP <- rhohat(mypatternPos,Z)
    
    
    ###################################
    Z <- distmap(mypatternStr)
    plot(Z,palette = 'Greens')
    myDist <- as.im(Z)
    myDist
    b = quantile(myDist, prob = (0:5)/5) # a way to quantify stroma interdistance!!!!!!
    bb = seq(from = 0, to = max(b), by = 25)
    bVS <- c()
    bVStotal <- c()
    myPr <- (0:5)/5
    plot(b,myPr)
    ####################################
    
  }
    
    write.table(myRhoP, file = paste(i,"_myRhoP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myRho, file = paste(i,"_myRhoN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    
    nn <- nncross(mypatternPos, mypatternStr)$dist
    mm <- nncross(mypatternNeg, mypatternStr)$dist
    write.table(nn, file = paste(i,"_myNNP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(mm, file = paste(i,"_myNNN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    
    ## starting PCF,J,G and L
    # RDF ###############################################################################################################################
    if (mypatternNeg$n > 50000){
      retain <- (runif(mypatternNeg$n) < 50000/mypatternNeg$n)
      thinNeg <- mypatternNeg[retain]}
    else{thinNeg <- mypatternNeg}
    
    
    XY <- superimpose(Pos=mypatternPos,Str=mypatternStr)
    XYN <- superimpose(Neg=thinNeg,Str=mypatternStr)
    PN <- superimpose(Pos=mypatternPos,Neg=thinNeg)
    
    mysim <- function(ppp1,ppp2){
      superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2)
    }
    
    mysimN <- function(ppp1,ppp2){
      superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2)
    }
    
    myenv <- envelope(XYN, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=4, simulate=expression(mysimN(PN,mypatternStr)))
    myenvp <- envelope(XY, fun = pcfcross, correction=c("Ripley"),nsim=39,r= rPCF, simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvp, file = paste(i,"_myenvRDFP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myenv, file = paste(i,"_myenvRDFN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    
    #Lcross ####################################################
    myenvL <- envelope(XYN, fun = Lcross, correction=c("Ripley"),nsim=4,r = rPCF, simulate=expression(mysimN(PN,mypatternStr)))
    myenvpL <- envelope(XY, fun = Lcross, correction=c("Ripley"),nsim=39, r= rPCF, simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvpL, file = paste(i,"_myenvLcrossP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myenvL, file = paste(i,"_myenvLcrossN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    #Jcross #########################################################
    XY <- superimpose(Pos=mypatternPos,Str=mypatternStr)
    XYN <- superimpose(Neg=mypatternNeg,Str=mypatternStr)
    PN <- superimpose(Pos=mypatternPos,Neg=mypatternNeg)
    
    mysim <- function(ppp1,ppp2){
      superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2)
    }
    
    mysimN <- function(ppp1,ppp2){
      superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2)
    }
    
    myenvJ <- envelope(XYN, fun = Jcross, correction=c("rs"),nsim=39,r =rJ, simulate=expression(mysimN(PN,mypatternStr)))
    myenvpJ <- envelope(XY, fun = Jcross, correction=c("rs"),nsim=39,r = rJ, simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvpJ, file = paste(i,"_myenvJcrossP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myenvJ, file = paste(i,"_myenvJcrossN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

    #Gcross #############################################################
    myenvG <- envelope(XYN, fun = Gcross, correction=c("rs"),nsim=39,r = rJ,  simulate=expression(mysimN(PN,mypatternStr)))
    myenvpG <- envelope(XY, fun = Gcross, correction=c("rs"),nsim=39, r = rJ,simulate=expression(mysim(PN,mypatternStr)))
    write.table(myenvpG, file = paste(i,"_myenvGcrossP",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    write.table(myenvG, file = paste(i,"_myenvGcrossN",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
    
  }}
