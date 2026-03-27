rm(list = ls())
library(sp)
library(rmapshaper)
library(ggplot2)
library(magrittr)
library(sf)
library(spatstat)
library(goftest)
library(readr)
library(kSamples)
library(maptools)
library(dplyr)
library(broom)
library(ggpubr)
library(ggplot2)

#setup and import data
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
myFile <- "A334.rds"#paste(readLines("fileName.txt"), collapse=" ")    
mynumber <- paste(readLines("number.txt"), collapse=" ") 
number <-as.numeric(mynumber)
results <- readRDS(myFile)

#Parameters
cell_diam <- 15 #13(Breast) #15(Lung) #microns
diameter <- TRUE #Whether to pixelate based on cell diameter 
scalefactor <- 0.5022 #micron per pixel

name <- names(results)
regions <- results[[name]]$regions
mydata <- read.csv(paste(getwd(),"/",number,"_Converted_coordinates_DIAMETER_",diameter,".csv",sep=""))

## get the individual large polygons as analysis windows
tumor_total_owin_gem <- as.data.frame(subset(regions,class_label=="tumor"))
polyList <- as.data.frame(tumor_total_owin_gem$geometry)

myArea <- c()
for (t in polyList$geometry){myArea <- c(myArea, st_area(t))}
myMax <- max(myArea)
mv_simpl <- st_simplify(subset(polyList$geometry, st_area(polyList$geometry) > myMax*0.25), preserveTopology = FALSE, 
                        dTolerance = 200)
plot(mv_simpl)

length(mv_simpl)
#No changes
nncrosss <- c()
nncrossr <- c()

qqVStotal <- rep(0, 8)
qqVRtotal <- rep(0, 8)
qqVSItotal <- rep(0, 8)
qqVRItotal <- rep(0, 8)
qqSNtotal <- rep(0, 8)# fraction of that type with respect to all cells within that interval
qqRNtotal <- rep(0, 8)

QStotal <- c()
QPtotal <- c()
QNtotal <- c()
bVStotal <- c()
numb <- 0

numPos <- 0 # number of sensitive cells in the sample
numNeg <- 0 # number of resistant cells in the sample
numEdgeTot <- 0 # number of tumor cells in stroma's vicinity
numStr <- 0 # number of stroma pixels in the sample
numEdgePos <- 0 # number of stroma in sensitive cells' vicinity
numEdgeNeg <- 0 # number of stroma in resisant cells' vicinity
numEdgeStr <- 0 #number of stroma edge pixels
numPosStrk1 <- 0 # number of sensitive tumor cells in stroma's immediate vicinity
numNegStrk1 <- 0 # number of resistant tumor cells in stroma's immediate vicinity

for (i in 1:length(mv_simpl)){
  
  all_owin <- as.owin(mv_simpl[i]*scalefactor)
  
  x <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Red')))]
  y <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Red')))]
  xy_inside<-inside.owin(x,y,all_owin)
  xx <-x[xy_inside]
  yy <- y[xy_inside]
  mypatternPos <- ppp(xx,yy,window=all_owin)#
  numPos <- numPos + length(mypatternPos$x)
  
  
  pdf(file = paste0(i, "_HeatContourSens.pdf"))
  K3 <- density(mypatternPos, sigma=150)
  plot(K3)
  contour(K3, add=TRUE)
  dev.off()
  
  
  xn <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Green')))]
  yn <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Green')))]
  xyn_inside<-inside.owin(xn,yn,all_owin)
  xxn <-xn[xyn_inside]
  yyn <- yn[xyn_inside]
  mypatternNeg <- ppp(xxn,yyn,window=all_owin)#
  numNeg <- numNeg + length(mypatternNeg$x)
  
  
  pdf(file = paste0(i, "_HeatContourRes.pdf"))
  K4 <- density(mypatternNeg, sigma=250)
  plot(K4)
  contour(K4, add=TRUE)
  dev.off()
  
  xs <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  ys <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  xys_inside<-inside.owin(xs,ys,all_owin)
  xxs <-xs[xys_inside]
  yys <- ys[xys_inside]
  mypatternStr <- ppp(xxs,yys,window=all_owin)#
  numStr <- numStr + length(mypatternStr$x)
  
  pdf(file = paste0(i, "_HeatContourSensStr.pdf"))
  K2 <- density(mypatternStr,sigma = 250)
  plot(K2)
  contour(K3, col = "green", add=TRUE)
  dev.off()
  
  pdf(file = paste0(i, "_HeatContourResStr.pdf"))
  K2 <- density(mypatternStr,sigma = 150)
  plot(K2)
  contour(K4, col = "green", add=TRUE)
  dev.off()

  pdf(file = paste0(i, "_HeatContourResStr.pdf"))
  #K2 <- density(mypatternPos,sigma = 250)
  plot(K3)
  contour(K2, col = "green", add=TRUE)
  dev.off()
  
  vx <- c()
  vy <- c()
  indexStr <- unique(nncross(mypatternNeg, mypatternStr, k =1)$which) 
  for (j in indexStr){
    vx <- c(vx, mypatternStr$x[j])
    vy <- c(vy, mypatternStr$y[j])}
  mypatternStrEdgeNeg <- ppp(vx, vy,window= all_owin )
  numEdgeNeg <- numEdgeNeg + length(mypatternStrEdgeNeg$x)
  
  vx <- c()
  vy <- c()
  indexStrp <- unique(nncross(mypatternPos, mypatternStr, k = 1)$which)
  for (h in indexStrp){
    vx <- c(vx, mypatternStr$x[h])
    vy <- c(vy, mypatternStr$y[h])}
  mypatternStrEdgePos <- ppp(vx, vy,window= all_owin )
  numEdgePos <- numEdgePos + length(mypatternStrEdgePos$x)
  
  xt <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Green', 'Red')))]
  yt <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Green', 'Red')))]
  xyt_inside<-inside.owin(xt,yt,all_owin)
  xxt <-xt[xyt_inside]
  yyt <- yt[xyt_inside]
  mypatternTot <- ppp(xxt,yyt,window=all_owin)#
  numEdgeTot <- numEdgeTot + length(mypatternTot$x)
  
  vx <- c()
  vy <- c()
  indexStrp <- unique(nncross(mypatternTot, mypatternStr, k = 2)$which)
  for (h in indexStrp){
    vx <- c(vx, mypatternStr$x[h])
    vy <- c(vy, mypatternStr$y[h])}
  mypatternStrEdgeTot <- ppp(vx, vy,window= all_owin )
  mypatternStrEdgeTot <- unique(mypatternStrEdgeTot)
  numEdgeStr <- numEdgeStr + length(mypatternStrEdgeTot$x)
  
  nnn <-nncross(mypatternPos, mypatternStrEdgeTot)$dist
  XPos_threshold_50_500 <- mypatternPos[nnn >= 50 & nnn <= 500]
  XPos_threshold_50_250 <- mypatternPos[nnn >= 50 & nnn <= 250]
  XPos_threshold_0_50 <- mypatternPos[nnn >= 0 & nnn <= 50]
  length(nnn)
  length(nnn[nnn>20])
  numPosStrk1 <- numPosStrk1 + length(nnn[nnn>20])
  mmm <- nncross(mypatternNeg, mypatternStrEdgeTot)$dist
  XNeg_threshold_50_500 <- mypatternNeg[mmm >= 50 & mmm <= 500]
  XNeg_threshold_50_250 <- mypatternNeg[mmm >= 50 & mmm <= 250]
  XNeg_threshold_0_50 <- mypatternNeg[mmm >= 0 & mmm <= 50]
  length(mmm)
  length(mmm[mmm>20])
  numNegStrk1 <- numNegStrk1 + length(mmm[mmm>20])
  
}

  mydataSum <- data.frame(numPos, 
                          numNeg, 
                          numEdgeTot, 
                          numStr, 
                          numEdgePos, 
                          numEdgeNeg, 
                          numEdgeStr, 
                          numPosStrk1, 
                          numNegStrk1)
  write.table(mydataSum, file = paste(number,"_SummaryStromaEdge",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
  
  
  Z <- distmap(mypatternStr)
  myDist <- as.im(Z)
  myTest <- auc(mypatternPos, myDist, high = TRUE)
  plot(myTestr <- roc(mypatternPos, myDist, high = TRUE))
  
  b = quantile(myDist, prob = (0:5)/5) # a way to quantify stroma interdistance!!!!!!
  bb = seq(from = 0, to = max(b), by = 25)
  bVS <- c()
  bVStotal <- c()
  myPr <- (0:5)/5
  
  
  plot(b,myPr, col = "black",  main = "Inter Stroma Distances",
       xlab = expression(bold(paste("Distance"," ", "(", mu,"m",")"))),
       ylab= "Cummulative Fraction of Distances ", ylim = c(0.000, 1.0),
       font.main=3, font.lab=2)
  
  

  
  # Build data frame from your vectors
  df <- data.frame(
    Distance = b,
    CDF = myPr
  )
  
  # Plot
  cairo_ps("CummStrFract.eps", width = 6, height = 4, family = "Arial")
  
  ggplot(df, aes(x = Distance, y = CDF)) +
    geom_point(shape = 21, fill = "white", color = "black", size = 2, stroke = 1) +
    labs(
      title = "Inter Stroma Distances",
      x = expression(bold("Distance")~"("~mu * "m"~")"),
      y = "Cummulative Fraction of Distances"
    ) +
    coord_cartesian(ylim = c(0, 1.0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 12)
    )
  
  ggsave(filename = "CummStrFract.eps", 
         device = cairo_ps, 
         dpi = 1200, 
         width = 6,
         height = 4, 
         family = "Arial")
  
  for (uu in 1:length(bb)){bVS <- c(bVS,bb[uu])}
  for (uuu in 1:length(bb)){bVStotal <- c(bVStotal,bb[uuu])}
  cuttest = cut(myDist, breaks = bb, label = bVStotal[1:(length(bb)-1)])#(1:(length(bb)-1)))
  testest = tess(image = cuttest)
  
  pdf(file = paste0(i, "_DistTess.pdf"))
  plot(testest)
  dev.off()
  
  qqS = quadratcount(mypatternPos, tess = testest)
  qqR = quadratcount(mypatternNeg, tess = testest)
  length(qqR)
  qqVS <- c() # fraction of the cells of that type within that interval
  qqVR <- c()
  qqVSI <- c()# intensity of that cell type within that interval
  qqVRI <- c()
  qqSN <- c()# fraction of that type with respect to all cells within that interval
  qqRN <- c()
  for (u in qqS){qqVS <- c(qqVS,u/mypatternPos$n)}
  for (v in qqR){qqVR <- c(qqVR,v/mypatternNeg$n)}
  for (uc in 1:length(qqS)){qqVSI <- c(qqVSI,qqS[uc]/area(testest[uc]))}
  for (uv in 1:length(qqR)){qqVRI <- c(qqVRI,qqR[uv]/area(testest[uv]))}
  for (uue in 1:length(qqS)){qqSN <- c(qqSN,qqS[uue]/(qqS[uue] + qqR[uue]))}
  for (uuu in 1:length(qqR)){qqRN <- c(qqRN,qqR[uuu]/(qqS[uuu] + qqR[uuu]))}
  
  pdf(file = paste0(i, "_Sensitive.pdf"))
  plot(bVS[1:length(qqVS)], qqVS)
  dev.off()
  pdf(file = paste0(i, "_Resistant.pdf"))
  plot(bVS[1:length(qqVR)], qqVR)
  dev.off()
  pdf(file = paste0(i, "_SensitiveNorm.pdf"))
  plot(bVS[1:length(qqVSI)], qqVSI)
  dev.off()
  pdf(file = paste0(i, "_ResistantNorm.pdf"))
  plot(bVS[1:length(qqVRI)], qqVRI)
  dev.off()
  
  pdf(file = paste0(i, "_ResistantNorm.pdf"))
  plot(bVS[1:8], qqSN[1:8], ylim = c(0.2,0.8), type = "b")
  #par(new=TRUE)
  points(bVS[1:8], qqRN[1:8], ylim= c(0.2,0.8), type = "b")
  dev.off()
  
  qqVStotal <- qqVStotal + qqVS[1:8]
  qqVRtotal <- qqVRtotal + qqVR[1:8]
  qqVSItotal <- qqVSItotal + qqVSI[1:8]
  qqVRItotal <- qqVRItotal + qqVRI[1:8]
  qqSNtotal <- qqSNtotal + qqSN[1:8]# fraction of that type with respect to all cells within that interval
  qqRNtotal <- qqRNtotal + qqRN[1:8]

  pdf(file= paste0(i,"_Plot.pdf"))
  plot(mypatternNeg, type = "p", cex = 0.00001, col = "red")
  plot(mypatternPos, type = "p", cex = 0.02, col = "green", add = TRUE)
  plot(mypatternStr, type = "p", cex = 0.02, col = "blue", add = TRUE)
  dev.off()

  
  # ss <- 5
  # QP <- quadratcount(mypatternPos, nx= ss, ny=ss)
  # QN <- quadratcount(mypatternNeg, nx= ss, ny=ss)
  # QS <- quadratcount(mypatternStr, nx= ss, ny=ss)
  # 
  # QS <- as.vector(t(QS))
  # QN <- as.vector(t(QN))
  # QP <- as.vector(t(QP))
  # my_data <- data.frame(QS, QP, QN)
  # 
  # for(jj in 1:length(QS)){
  #   QStotal <- c(QStotal, QS[jj])}
  # 
  # for(jj in 1:length(QP)){
  #   QPtotal <- c(QPtotal, QP[jj]/(QP[jj] + QN[jj]))}
  # 
  # for(jj in 1:length(QN)){
  #   QNtotal <- c(QNtotal, QN[jj]/(QP[jj] + QN[jj]))}
  # 
  # library("ggpubr")
  # pdf(file = paste0(i, "StromaSensitiveCorr.pdf"))
  # ggscatter(my_data, x = "QS", y = "QP", 
  #           add = "reg.line", conf.int = TRUE, 
  #           cor.coef = TRUE, cor.method = "spearman",
  #           xlab = "Number of Stroma Pixels", ylab = "Number Sensitive Cells")
  # dev.off()
  # 
  # pdf(file = paste0(i, "StromaResistantCorr.pdf"))
  # ggscatter(my_data, x = "QS", y = "QN", 
  #           add = "reg.line", conf.int = TRUE, 
  #           cor.coef = TRUE, cor.method = "spearman",
  #           xlab = "Number of Stroma Pixels", ylab = "Number Resistant Cells")
  # dev.off()
  # 
  # corSP <- cor(QS, QP, method = c("pearson", "kendall", "spearman"))
  # corSPtest <- cor.test(QS, QP, method=c("pearson", "kendall", "spearman"))
  
  # check to not analize too small off polygons
 
    
    nn <- nncross(mypatternPos, mypatternStr)$dist
    mm <- nncross(mypatternNeg, mypatternStr)$dist
    for(jj in 1:length(nn)){
      nncrosss <- c(nncrosss, nn[jj])}
    for(hh in 1:length(mm)){
      nncrossr <- c(nncrossr, mm[hh])}
    

nn_0_50SR <- nncross(XPos_threshold_0_50, XNeg_threshold_0_50)$dist
plot(density(nn_0_50SR))
nn_0_50SStr <- nncross(XPos_threshold_0_50, mypatternStr)$dist
plot(density(nn_0_50SStr))
nn_0_50RStr <- nncross(XNeg_threshold_0_50, mypatternStr)$dist
plot(density(nn_0_50RStr))


plot(density(nn_0_50SR), xlim = c(0, 250), ylim = c(0,0.035), 
     col = "black", lwd = 3, 
     xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Density", font.lab=2, label = "RS")
lines(density(nn_0_50SStr),col = "darkgreen", lwd = 3,label = "SStr")
lines(density(nn_0_50RStr),col = "red", lwd = 3,label = "RStr")
legend()


#######ggplot
# Compute density objects
dens_RS   <- density(nn_0_50SR)
dens_SStr <- density(nn_0_50SStr)
dens_RStr <- density(nn_0_50RStr)

# Combine into one long-format data frame
df <- data.frame(
  x = c(dens_RS$x, dens_SStr$x, dens_RStr$x),
  y = c(dens_RS$y, dens_SStr$y, dens_RStr$y),
  group = rep(c("RS", "SStr", "RStr"), each = length(dens_RS$x))
)


cairo_ps("334_Dist0_50.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("RS" = "black", "SStr" = "darkgreen", "RStr" = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.035)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )

ggsave(filename = "334_Dist0_50.eps.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
################

nn_50_250SR <- nncross(XPos_threshold_50_250, XNeg_threshold_50_250)$dist
plot(density(nn_50_250SR))
nn_50_250SStr <- nncross(XPos_threshold_50_250, mypatternStr)$dist
plot(density(nn_50_250SStr))
nn_50_250RStr <- nncross(XNeg_threshold_50_250, mypatternStr)$dist
plot(density(nn_50_250RStr))

#######ggplot
# Compute density objects
dens_RS   <- density(nn_50_250SR)
dens_SStr <- density(nn_50_250SStr)
dens_RStr <- density(nn_50_250RStr)

# Combine into one long-format data frame
df <- data.frame(
  x = c(dens_RS$x, dens_SStr$x, dens_RStr$x),
  y = c(dens_RS$y, dens_SStr$y, dens_RStr$y),
  group = rep(c("RS", "SStr", "RStr"), each = length(dens_RS$x))
)


cairo_ps("334_Dist50_250.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("RS" = "black", "SStr" = "darkgreen", "RStr" = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.035)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )

ggsave(filename = "334_Dist50_250.eps.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
################


nn_50_500SR <- nncross(XPos_threshold_50_500, XNeg_threshold_50_500)$dist
plot(density(nn_50_500SR))
nn_50_500SStr <- nncross(XPos_threshold_50_500, mypatternStr)$dist
plot(density(nn_50_500SStr))
nn_50_500RStr <- nncross(XNeg_threshold_50_500, mypatternStr)$dist
plot(density(nn_50_500RStr))

#######ggplot
# Compute density objects
dens_RS   <- density(nn_50_500SR)
dens_SStr <- density(nn_50_500SStr)
dens_RStr <- density(nn_50_500RStr)

# Combine into one long-format data frame
df <- data.frame(
  x = c(dens_RS$x, dens_SStr$x, dens_RStr$x),
  y = c(dens_RS$y, dens_SStr$y, dens_RStr$y),
  group = rep(c("RS", "SStr", "RStr"), each = length(dens_RS$x))
)


cairo_ps("334_Dist50_500.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("RS" = "black", "SStr" = "darkgreen", "RStr" = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.035)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )

ggsave(filename = "334_Dist50_500.eps.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
################

write.table(nncrosss, file = paste(number,"_myNNS",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(nncrossr, file = paste(number,"_myNNR",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

write.table(QPtotal, file = paste(number,"_myQPtotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(QStotal, file = paste(number,"_myQStotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(QNtotal, file = paste(number,"_myQNtotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(bVStotal, file = paste(number,"_mybVStotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

for (tt  in 1:8){
  qqVStotal[tt] <- qqVStotal[tt]/length(mv_simpl)
  qqVRtotal[tt] <- qqVRtotal[tt]/length(mv_simpl)
  qqVSItotal[tt] <- qqVSItotal[tt]/length(mv_simpl)
  qqVRItotal[tt] <- qqVRItotal[tt]/length(mv_simpl)
  qqSNtotal[tt] <- qqSNtotal[tt]/length(mv_simpl)# fraction of that type with respect to all cells within that interval
  qqRNtotal[tt] <- qqRNtotal[tt]/length(mv_simpl)}


write.table(qqVStotal, file = paste(number,"_myqqVStotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(qqVRtotal, file = paste(number,"_myqqVRtotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(qqVSItotal, file = paste(number,"_myqqVSItotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(qqVRItotal, file = paste(number,"_myqqVRItotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(qqSNtotal, file = paste(number,"_myqqSNtotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")
write.table(qqRNtotal, file = paste(number,"_myqqRNtotal",".csv",sep=""),row.names=FALSE,col.names=TRUE,sep=",")

# my_data_total <- data.frame(QStotal, QPtotal, QNtotal)
# pdf(file = paste0(i, "StrSensCorrtotal.pdf"))
# ggscatter(my_data_total, x = "QStotal", y = "QPtotal", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman",
#           xlab = "Number of Stroma Pixels", ylab = "Number Sensitive Cells")
# dev.off()
# 
# 
# 
# pdf(file = paste0(i, "StrResistCorrtotal.pdf"))
# ggscatter(my_data_total, x = "QStotal", y = "QNtotal", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "spearman",
#           xlab = "Number of Stroma Pixels", ylab = "Number Resistant Cells")
# dev.off()



nn <- density(nncrosss)
mm <- density(nncrossr)
xmaxs <- which.max(density(nncrosss)$y)
density(nncrosss)$x[xmaxs]
xmaxr <- which.max(density(nncrossr)$y)
density(nncrossr)$x[xmaxr]


pdf(file=paste0(number,"_DistToStromaWholeSample.pdf"),width = 8, height = 8)
plot(density(nncrosss), col = "darkgreen", lwd = 3, xlim = c(-0, 150),
     xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Density", font.lab=2)
lines(density(nncrossr), col = "red", lwd = 3)
polygon(nn, density = -1, col = rgb(0.4,0.6,0.2, 0.5))  
polygon(mm, density = -1, col = rgb(0.9, 0, 0.1, 0.4))
abline(v = density(nncrosss)$x[xmaxs], col = "darkgreen", lty = 1, lwd = 2)
abline(v = density(nncrossr)$x[xmaxr], col = "red", lty = 1, lwd = 2)
dev.off()



################### ggplot
library(ggplot2)

# Compute densities
dens_s <- density(nncrosss)
dens_r <- density(nncrossr)

# Create data frames for plotting
df_s <- data.frame(x = dens_s$x, y = dens_s$y, group = "Sensitive")
df_r <- data.frame(x = dens_r$x, y = dens_r$y, group = "Resistant")
df_density <- rbind(df_s, df_r)

# Create vertical lines at peaks
peak_lines <- data.frame(
  x = c(dens_s$x[xmaxs], dens_r$x[xmaxr]),
  group = c("Sensitive", "Resistant"),
  color = c("darkgreen", "red")
)

# Assuming nn and mm are data frames or matrices with x and y for the shaded areas
# If they are just numeric vectors, wrap them accordingly
#df_nn <- as.data.frame(nn)
#colnames(df_nn) <- c("x", "y")

#df_mm <- as.data.frame(mm)
#colnames(df_mm) <- c("x", "y")

# Plot
cairo_ps("334_DistToStromaWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot() +
  # Shaded polygons
  #geom_polygon(data = df_nn, aes(x = x, y = y), fill = rgb(0.4, 0.6, 0.2, 0.5)) +
  #geom_polygon(data = df_mm, aes(x = x, y = y), fill = rgb(0.9, 0, 0.1, 0.4)) +
  
  # Density lines
  geom_line(data = df_density, aes(x = x, y = y, color = group), size = 1.3) +
  
  # Vertical lines at density peaks
  #geom_vline(data = peak_lines, aes(xintercept = x, color = group), linetype = "solid", linewidth = 1) +
  
  scale_color_manual(values = c("Sensitive" = "darkgreen", "Resistant" = "red")) +
  
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density"
  ) +
  coord_cartesian(xlim = c(0, 150)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13, vjust= 5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20),
    axis.text = element_text(size = 12),
    legend.position = "top"
  )

ggsave(filename = "334_DistToStromaWholeSample.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
######################

setEPS()
postscript(file=paste0(number,"_DistToStromaWholeSample2.eps"),width = 8, height = 8)
plot(density(nncrosss), col = "darkgreen", lwd = 3, xlim = c(-0, 150),
     xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Density", font.lab=2)
lines(density(nncrossr), col = "red", lwd = 3)
polygon(nn, density = -1, col = rgb(0.4,0.6,0.2, 0.7))  
polygon(mm, density = -1, col = rgb(0.9, 0, 0.1, 0.4))
#abline(v = density(nncrosss)$x[xmaxs], col = "darkgreen", lty = 1, lwd = 2)
#abline(v = density(nncrossr)$x[xmaxr], col = "red", lty = 1, lwd = 2)
dev.off()

f <- ecdf(nncrosss)
ff <- ecdf(nncrossr)
myKS <- ks.test(nncrosss, ff)
pdf(file=paste0(number,"_CumKSWholeSample.pdf"),width = 8, height = 8)
plot(f, col = "darkgreen", main = "NN Cumulative Density Function",xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     xlim = c(0,150),
     ylab= "Sensitive/Resistant Cumulative Density",
     font.main=3, font.lab=2, lty = 1, lwd = 3)
lines(ff, col = "red", lwd = 3)
text(120, 0.2, "KS  = 0.296")
dev.off()


##########ggplot
# Convert step function to data frame
df_f <- data.frame(x = knots(f), y = f(knots(f)), group = "Sensitive")
df_ff <- data.frame(x = knots(ff), y = ff(knots(ff)), group = "Resistant")
df_all <- rbind(df_f, df_ff)

cairo_ps("334_CumKSWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot(df_all, aes(x = x, y = y, color = group)) +
  geom_step(size = 1.3) +  # Use step to mimic base R's ecdf style
  scale_color_manual(values = c("Sensitive" = "darkgreen", "Resistant" = "red")) +
  labs(
    title = "NN Cumulative Density Function",
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Cumulative Density"
  ) +
  coord_cartesian(xlim = c(0, 150)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13, vjust= 5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20),
    axis.text = element_text(size = 12),
    legend.position = "top"
  )
ggsave(filename = "334_CumKSWholeSample.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

postscript(file=paste0(number,"_CumKSWholeSample.eps"),width = 8, height = 8)
plot(f, col = "darkgreen", main = "NN Cumulative Density Function",xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     xlim = c(0,150),
     ylab= "Sensitive/Resistant Cumulative Density",
     font.main=3, font.lab=2, lty = 1, lwd = 3)
lines(ff, col = "red", lwd = 3)
text(120, 0.2, "KS  = 0.296")
dev.off()

#plot(bVStotal, qqVStotal, col = "red", type = "b", lwd = 3)
#points(bVStotal, qqVRtotal, col = "darkgreen", type = "b", lwd = 3)
pdf(file=paste0(number,"_SensResistanceDens.pdf"),width = 8, height = 8)
plot(bVS[1:8], qqVRItotal[1:8], col = "red", type = "b", lwd = 3, main = "Fractions of Resistance and Sensitive",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Densities", ylim = c(0.00001, 0.00055),
     font.main=3, font.lab=2)
points(bVS[1:8], qqVSItotal[1:8], col = "darkgreen", type = "b", lwd = 3)
dev.off()

pdf(file=paste0(number,"_SensResistanceNumb.pdf"),width = 8, height = 8)
plot(bVS[1:8], qqVRtotal[1:8], col = "red", type = "b", lwd = 3, main = "Fractions of Resistance and Sensitive",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Fractions", ylim = c(0.000, 0.45),
     font.main=3, font.lab=2)
points(bVS[1:8], qqVStotal[1:8], col = "darkgreen", type = "b", lwd = 3)
dev.off()


pdf(file=paste0(number,"_SensResistanceFract.pdf"),width = 8, height = 8)
plot(bVS[1:8], qqRNtotal[1:8], col = "red", type = "b", lwd = 3, main = "Fractions of Resistance and Sensitive",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Fractions", ylim = c(0.000, 1.0),
     font.main=3, font.lab=2)
points(bVS[1:8], qqSNtotal[1:8], col = "darkgreen", type = "b", lwd = 3)
dev.off()

Z <- distmap(mypatternStr)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,Z)
plot(myRho[1:100])
plot(myRhoP[1:100])
rLim <- 210


postscript(file=paste0(number, "_rhohatcombinedtwoAxesWSample.eps"))
par(mar = c(5, 4, 4, 4) + 0.3)              
plot(myRhoP$Z[1:rLim], myRhoP$rho[1:rLim],  main="Rho hat BrdU+ and BrdU-",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     #ylim= c(min(myRhoP$lo[1:rLim]),max(myRhoP$hi[1:rLim])),
     type='l', col= "darkgreen", lwd = 3,
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "darkgreen")     
lines(myRhoP$Z[1:rLim], myRhoP$hi[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rLim], myRhoP$lo[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rLim], rev(myRhoP$Z[1:rLim])), c(myRhoP$lo[1:rLim], rev(myRhoP$hi[1:rLim]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
par(new = TRUE)  
plot(myRho$Z[1:rLim], myRho$rho[1:rLim], type='l', col= "red", lwd = 3, 
     ylim= c(min(myRho$lo),max(myRho$hi)),   bty = "n",       
     axes = FALSE, xlab = "", ylab = "")
lines(myRho$Z[1:rLim], myRho$hi[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRho$Z[1:rLim], myRho$lo[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRho$Z[1:rLim], rev(myRho$Z[1:rLim])), c(myRho$lo[1:rLim], rev(myRho$hi[1:rLim]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRho$ave[20], lty = 2, lwd = 2, col = rgb(0.9,0.2,0))
axis(side = 4)#, at = pretty(range(myRho$hi[1:rLim]))) 
mtext("Resistant Tumor Cells RSF",  side = 4, line=2.5, cex=1.3, font=2, col ="red" )
dev.off()

scaling_factor <- max(myRhoP$hi[1:rLim]) / max(myRho$hi[1:rLim])

# Set the limits for the second axis based on the scaled values of y2
y2_scaled_range <- c(0, max(myRho$hi[1:rLim]) * scaling_factor)
pdf(file=paste0(number, "_rhohatcombinedtwoAxesWSample.pdf"))
par(mar = c(5, 4, 4, 4) + 0.3)              
plot(myRhoP$Z[1:rLim], myRhoP$rho[1:rLim],  main="Rho hat BrdU+ and BrdU-",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     #ylim= c(min(myRhoP$lo[1:rLim]),max(myRhoP$hi[1:rLim])),
     type='l', col= "darkgreen", lwd = 3,
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "darkgreen")     
lines(myRhoP$Z[1:rLim], myRhoP$hi[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rLim], myRhoP$lo[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rLim], rev(myRhoP$Z[1:rLim])), c(myRhoP$lo[1:rLim], rev(myRhoP$hi[1:rLim]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
par(new = TRUE)  
plot(myRho$Z[1:rLim], myRho$rho[1:rLim]* scaling_factor, type='l', col= "red", lwd = 3, 
     ylim= c(min(myRho$lo),max(myRho$hi)),   bty = "n",       
     axes = FALSE, xlab = "", ylab = "")
lines(myRho$Z[1:rLim], myRho$hi[1:rLim]* scaling_factor, col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRho$Z[1:rLim], myRho$lo[1:rLim]* scaling_factor, col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRho$Z[1:rLim], rev(myRho$Z[1:rLim])), 
        c(myRho$lo[1:rLim]* scaling_factor, rev(myRho$hi[1:rLim]* scaling_factor)),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRho$ave[20], lty = 2, lwd = 2, col = rgb(0.9,0.2,0))
axis(side=4)#, at=seq(0, max(myRho$hi[1:rLim]), length.out=5), labels=seq(0, max(myRho$hi[1:rLim]), length.out=5)* scaling_factor, las=1)
#axis(side = 4)#, at = pretty(range(0, myRho$hi[1:rLim]))) 
mtext("Resistant Tumor Cells RSF",  side = 4, line=2.5, cex=1.3, font=2, col ="red" )
dev.off()

rRho <-215
pdf(file=paste0(1, "_rhohatcombinedtwoAxesSample.pdf"))
par(mar = c(5, 4, 4, 4) + 0.3)              
plot(myRhoP$Z[1:rRho], myRhoP$rho[1:rRho],  main="Rho hat Sensitive and Resistant Tumor Cells",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     ylim= c(min(myRhoP$lo[1:rRho]),max(myRhoP$hi[1:rRho])),
     type='l', col= "red", lwd = 3,
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "red")     
lines(myRhoP$Z[1:rRho], myRhoP$hi[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rRho], myRhoP$lo[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rRho], rev(myRhoP$Z[1:rRho])), c(myRhoP$lo[1:rRho], rev(myRhoP$hi[1:rRho]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
par(new = TRUE)  
plot(myRho$Z[1:rRho], myRho$rho[1:rRho], type='l', col= "darkgreen", lwd = 3, 
     ylim= c(min(myRho$lo[1:(rRho+10)]),max(myRho$hi[1:(rRho+10)])),   bty = "n",       
     axes = FALSE, xlab = "", ylab = "")
lines(myRho$Z[1:rRho], myRho$hi[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRho$Z[1:rRho], myRho$lo[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRho$Z[1:rRho], rev(myRho$Z[1:rRho])), c(myRho$lo[1:rRho], rev(myRho$hi[1:rRho]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRho$ave[20], lty = 2, lwd = 2, col = rgb(0.9,0.2,0))
axis(side = 4)#, at = pretty(range(myRho$hi[1:150])))#at = pretty(range(myRho$hi[1:250]))) 
mtext("Resistant Tumor Cells RSF",  side = 4, line=2.5, cex=1.3, font=2, col ="darkgreen" )
dev.off()

postscript(file=paste0(number, "_rhohatcombinedtwoAxesWSample.eps"))
par(mar = c(5, 4, 4, 4) + 0.3)              
plot(myRhoP$Z[1:rRho], myRhoP$rho[1:rRho],  main="Rho hat Sensitive and Resistant Tumor Cells",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     ylim= c(min(myRhoP$lo[1:rRho]),max(myRhoP$hi[1:rRho])),
     type='l', col= "red", lwd = 3,
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "red")     
lines(myRhoP$Z[1:rRho], myRhoP$hi[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rRho], myRhoP$lo[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rRho], rev(myRhoP$Z[1:rRho])), c(myRhoP$lo[1:rRho], rev(myRhoP$hi[1:rRho]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
par(new = TRUE)  
plot(myRho$Z[1:rRho], myRho$rho[1:rRho], type='l', col= "darkgreen", lwd = 3, 
     ylim= c(min(myRho$lo[1:(rRho+10)]),max(myRho$hi[1:(rRho+10)])),   bty = "n",       
     axes = FALSE, xlab = "", ylab = "")
lines(myRho$Z[1:rRho], myRho$hi[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRho$Z[1:rRho], myRho$lo[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRho$Z[1:rRho], rev(myRho$Z[1:rRho])), c(myRho$lo[1:rRho], rev(myRho$hi[1:rRho]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRho$ave[20], lty = 2, lwd = 2, col = rgb(0.9,0.2,0))
axis(side = 4)#, at = pretty(range(myRho$hi[1:150])))#at = pretty(range(myRho$hi[1:250]))) 
mtext("Resistant Tumor Cells RSF",  side = 4, line=2.5, cex=1.3, font=2, col ="darkgreen" )
dev.off()

################ggplot
install.packages("scales")  # if you haven't already
library(scales)
# Get min and max for y-values
Z <- distmap(mypatternStr)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,Z)
plot(myRho[1:100])
plot(myRhoP[1:100])
rLim <- 210
rRho <-215
range_sensitive <- range(myRhoP$rho[1:rRho])
range_resistant <- range(myRho$rho[1:rRho])

# Rescale function: resistant -> sensitive scale
scale_factor <- diff(range_sensitive) / diff(range_resistant)
offset <- range_sensitive[1] - range_resistant[1] * scale_factor

# Apply to resistant data
scaled_rho <- myRho$rho[1:rRho] * scale_factor + offset
scaled_hi <- myRho$hi[1:rRho] * scale_factor + offset
scaled_lo <- myRho$lo[1:rRho] * scale_factor + offset
scaled_ave <- myRho$ave[20] * scale_factor + offset

Z_vals <- myRhoP$Z[1:rRho]

df <- data.frame(
  Z = c(Z_vals, Z_vals),
  rho = c(myRhoP$rho[1:rRho], scaled_rho),
  lo = c(myRhoP$lo[1:rRho], scaled_lo),
  hi = c(myRhoP$hi[1:rRho], scaled_hi),
  group = rep(c("Sensitive", "Resistant"), each = rRho)
)

ref_lines <- data.frame(
  yintercept = c(myRhoP$ave[20], scaled_ave),
  group = c("Sensitive", "Resistant")
)

cairo_ps("334_rhohatCombinedtwoAxesWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot(df, aes(x = Z, y = rho, group = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  geom_line(aes(color = group), size = 1.2) +
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  geom_hline(data = ref_lines, aes(yintercept = yintercept, color = group),
             linetype = "dashed", size = 1) +
  
  scale_color_manual(values = c("Sensitive" = "darkgreen", "Resistant" = "red")) +
  
  scale_y_continuous(
    name = "Sensitive Tumor Cells RSF",
    labels = scales::scientific,
    sec.axis = sec_axis(
      ~ (. - offset) / scale_factor,
      name = "Resistant Tumor Cells RSF",
      labels = scales::scientific
    )
  ) +
  
  labs(
    title = "Rho hat Sensitive and Resistant Tumor Cells",
    x = expression(bold(paste("Distance from Stroma (", mu, "m)")))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove full plot border
    axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
    axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
    axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
    axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
    
    
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y.left = element_text(face = "bold", size = 13, color = "darkgreen"),
    axis.title.y.right = element_text(face = "bold", size = 13, color = "red"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
ggsave(filename = "334_rhohatCombinedtwoAxesWholeSample.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")


############ggplot
Z <- distmap(mypatternNeg)
myRhoP <- rhohat(XPos_threshold_0_50,Z)
plot(myRhoP$Z, myRhoP$rho, main="Rho hat Sensitive Tumor Cells",
     xlab = expression(bold(paste("Distance from Resistant Tumor Cells"," ", "(", mu,"m",")"))),
     type='l', col= "black", lwd = 3, xlim = c(0,155),
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "black")   
lines(myRhoP$Z[1:rRho], myRhoP$hi[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rRho], myRhoP$lo[1:rRho], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rRho], rev(myRhoP$Z[1:rRho])), c(myRhoP$lo[1:rRho], rev(myRhoP$hi[1:rRho]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = "black")

df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
cairo_ps("334_rhohatSRWholeSample0_50.eps", width = 8, height = 6, family = "Arial")
ggplot(df[1:rRho, ], aes(x = Z)) +
  # Shaded confidence interval
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  
  # Main RSF line
  geom_line(aes(y = rho), color = "black", size = 1.3) +
  
  # Confidence bounds
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  
  # Horizontal reference line
  geom_hline(yintercept = myRhoP$ave[20], linetype = "dashed", size = 1.2, color = "black") +
  
  labs(
    title = "Rho hat Sensitive Tumor Cells",
    x = expression(bold(paste("Distance from Resistant Tumor Cells (", mu, "m)"))),
    y = "Sensitive tumor Cells RSF"
  ) +
  coord_cartesian(xlim = c(0, 155)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        
    plot.title = element_text(face = "italic", size = 14, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12)
  )

ggsave(filename = "334_rhohatSRWholeSample0_50.eps",
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

#################


#######ggplot
Z <- distmap(mypatternNeg)
myRhoP <- rhohat(XPos_threshold_50_250,Z)
df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
cairo_ps("334_rhohatSRWholeSample50_250.eps", width = 8, height = 6, family = "Arial")
ggplot(df[1:rRho, ], aes(x = Z)) +
  # Shaded confidence interval
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  
  # Main RSF line
  geom_line(aes(y = rho), color = "black", size = 1.3) +
  
  # Confidence bounds
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  
  # Horizontal reference line
  geom_hline(yintercept = myRhoP$ave[20], linetype = "dashed", size = 1.2, color = "black") +
  
  labs(
    title = "Rho hat Sensitive Tumor Cells",
    x = expression(bold(paste("Distance from Resistant Tumor Cells (", mu, "m)"))),
    y = "Sensitive tumor Cells RSF"
  ) +
  coord_cartesian(xlim = c(0, 155)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        
        plot.title = element_text(face = "italic", size = 14, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12)
  )

ggsave(filename = "334_rhohatSRWholeSample50_250.eps",
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

#################
#######ggplot
Z <- distmap(mypatternNeg)
myRhoP <- rhohat(XPos_threshold_50_500,Z)
df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
cairo_ps("334_rhohatSRWholeSample50_500.eps", width = 8, height = 6, family = "Arial")
ggplot(df[1:rRho, ], aes(x = Z)) +
  # Shaded confidence interval
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  
  # Main RSF line
  geom_line(aes(y = rho), color = "black", size = 1.3) +
  
  # Confidence bounds
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  
  # Horizontal reference line
  geom_hline(yintercept = myRhoP$ave[20], linetype = "dashed", size = 1.2, color = "black") +
  
  labs(
    title = "Rho hat Sensitive Tumor Cells",
    x = expression(bold(paste("Distance from Resistant Tumor Cells (", mu, "m)"))),
    y = "Sensitive tumor Cells RSF"
  ) +
  coord_cartesian(xlim = c(0, 155)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        
        plot.title = element_text(face = "italic", size = 14, hjust = 0.5),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12)
  )

ggsave(filename = "334_rhohatSRWholeSample50_500.eps",
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

#################



## starting PCF,J,G and L
# plot the pcf or rdf
XY <- superimpose(Pos=mypatternPos,Str=mypatternStr)
XYN <- superimpose(Neg=mypatternNeg,Str=mypatternStr)
PN <- superimpose(Pos=mypatternPos,Neg=mypatternNeg)
PN_0_50 <- superimpose(Pos=XPos_threshold_0_50,Neg=XNeg_threshold_0_50)
PN_50_250 <- superimpose(Pos=XPos_threshold_50_250,Neg=XNeg_threshold_50_250)
PN_50_500 <- superimpose(Pos=XPos_threshold_50_500,Neg=XNeg_threshold_50_500)
mysim <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2)
}

mysimN <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2)
}

rPCF <-seq(from = 0, to = 450, by = 0.5)
myCrossEnv <- envelope(PN, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=39, simulate=expression(mysimN(PN,mypatternStr)))
plot (myCrossEnv$r[19:550],(myCrossEnv$obs/myCrossEnv$mmean)[19:550])

pdf(file=paste0(number, "_PCFSenResWholeSample.pdf"))
plot(myCrossEnv$r[19:650],(myCrossEnv$obs/myCrossEnv$mmean)[19:650],type = 'l',main = "PCF Sensitive and Resistant",
     xlab = expression(bold(paste("Cell-Cell Distance"," ", "(", mu,"m",")"))), 
     cex.lab = 1.3, font.lab = 4,  ylim = c(0.0, 10.5), xlim = c(10, 320), col = "black", lwd = 4,
     ylab =""# expression(bold("g"["BrdU+,Str"])),
)
title(ylab=expression(bold("g"["Sensitive,Resistant"])), line=2.5, cex.lab=1.4)
abline(h = 1)
dev.off()


myCrossEnv_0_50 <- envelope(PN_0_50, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=39, simulate=expression(mysimN(PN_0_50,mypatternStr)))
plot(myCrossEnv_0_50)

myCrossEnv_50_250 <- envelope(PN_50_250, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=39, simulate=expression(mysimN(PN_50_250,mypatternStr)))
plot(myCrossEnv_50_250)

myCrossEnv_50_500 <- envelope(PN_50_500, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=39, simulate=expression(mysimN(PN_50_500,mypatternStr)))
plot(myCrossEnv_50_500)


postscript(file=paste0(number, "_PCFSenResWholeSample.eps"))
plot(myCrossEnv$r[19:650],(myCrossEnv$obs/myCrossEnv$mmean)[19:650],type = 'l',main = "PCF Sensitive and Resistant",
     xlab = expression(bold(paste("Cell-Cell Distance"," ", "(", mu,"m",")"))), 
     cex.lab = 1.3, font.lab = 4,  ylim = c(0.0, 10.5), xlim = c(10, 320), col = "black", lwd = 4,
     ylab =""# expression(bold("g"["BrdU+,Str"])),
)
title(ylab=expression(bold("g"["Sensitive,Resistant"])), line=2.5, cex.lab=1.4)
abline(h = 1)
dev.off()

grDevices::cairo_ps()
# Prepare the data
df_pcf <- data.frame(
  r = myCrossEnv$r[0:650],
  g = (myCrossEnv$obs / myCrossEnv$mmean)[0:650]
)
cairo_ps("RDSSenRes.eps", width = 6, height = 4, family = "Arial")

# Plot using ggplot2
ggplot(df_pcf, aes(x = r, y = g)) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 320)) +
  scale_y_continuous(limits = c(0.0, 10.5)) +
  labs(
    title = "PCF Sensitive and Resistant",
    x = expression(bold(paste("Cell-Cell Distance", " ", "(", mu, "m", ")"))),
    y = expression(bold("g"["Sensitive,Resistant"]))
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11)
  )

dev.off()


df_pcf <- data.frame(
  r = myCrossEnv_0_50$r[0:650],
  g = (myCrossEnv_0_50$obs / myCrossEnv_0_50$mmean)[0:650]
)
cairo_ps("RDSSenRes_0_50.eps", width = 6, height = 4, family = "Arial")
# Plot using ggplot2
ggplot(df_pcf, aes(x = r, y = g)) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 320)) +
  scale_y_continuous(limits = c(0.0, 10.5)) +
  labs(
    title = "PCF Sensitive and Resistant",
    x = expression(bold(paste("Cell-Cell Distance", " ", "(", mu, "m", ")"))),
    y = expression(bold("g"["Sensitive,Resistant"]))
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 11)
  )
dev.off()


df_pcf <- data.frame(
  r = myCrossEnv_50_250$r[0:650],
  g = (myCrossEnv_50_250$obs / myCrossEnv_50_250$mmean)[0:650]
)
cairo_ps("RDSSenRes_50_250.eps", width = 6, height = 4, family = "Arial")
# Plot using ggplot2
ggplot(df_pcf, aes(x = r, y = g)) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 320)) +
  scale_y_continuous(limits = c(0.0, 10.5)) +
  labs(
    title = "PCF Sensitive and Resistant",
    x = expression(bold(paste("Cell-Cell Distance", " ", "(", mu, "m", ")"))),
    y = expression(bold("g"["Sensitive,Resistant"]))
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 11)
  )
dev.off()


df_pcf <- data.frame(
  r = myCrossEnv_50_500$r[0:650],
  g = (myCrossEnv_50_500$obs / myCrossEnv_50_500$mmean)[0:650]
)
cairo_ps("RDSSenRes_50_500.eps", width = 6, height = 4, family = "Arial")
# Plot using ggplot2
ggplot(df_pcf, aes(x = r, y = g)) +
  geom_line(color = "black", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 320)) +
  scale_y_continuous(limits = c(0.0, 10.5)) +
  labs(
    title = "PCF Sensitive and Resistant",
    x = expression(bold(paste("Cell-Cell Distance", " ", "(", mu, "m", ")"))),
    y = expression(bold("g"["Sensitive,Resistant"]))
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),  # Remove full plot border
        axis.line.x = element_line(color = "black", size = 1),  # Keep bottom axis
        axis.line.y = element_line(color = "black", size = 1),  # Keep left axis
        axis.ticks.y = element_line(color = "black", size = 1),  # Keep Y-axis ticks
        axis.ticks.x = element_line(color = "black", size = 1),#element_blank(),  # Remove X-axis ticks
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.title.y = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 11)
  )
dev.off()



myenv <- envelope(XYN, fun = pcfcross, correction=c("Ripley"),r = rPCF, nsim=39, simulate=expression(mysimN(PN,mypatternStr)))
myenvp <- envelope(XY, fun = pcfcross, correction=c("Ripley"),nsim=39,r= rPCF, simulate=expression(mysim(PN,mypatternStr)))

plot(myenvp)
plot(myenv)

pdf(file=paste0(number, "_PCFStrSenResWholeSample.pdf"))
plot(myenvp$r[19:650],(myenvp$obs/myenvp$mmean)[19:650],type = 'l',main = "PCF Sensitive and Resistant",
     xlab = expression(bold(paste("Stroma-Cell Distance"," ", "(", mu,"m",")"))), 
     cex.lab = 1.3, font.lab = 4,  ylim = c(0.0, 1.5), xlim = c(10, 320), col = "darkgreen", lwd = 3,
     ylab =""# expression(bold("g"["BrdU+,Str"])),
)
title(ylab=expression(bold("g"["Stroma,Sensitive/Resistant"])), line=2.5, cex.lab=1.4)
lines(myenv$r[19:650],(myenv$obs/myenv$mmean)[19:650],col = "red", lty = 1, lwd = 3)
# lines(x,(averObs)[15:PCF_r_length], col = rgb(0.3,0.2,1), lty = 1, lwd = 3)
# lines(x, y1, col = rgb(0.5,0.5,0.5), lty = 2, lwd = 1)
lines(myenvp$r[19:650],(myenvp$mmean/myenvp$mmean)[19:650], col = rgb(0.5,0.5,0.5), lty = 2, lwd = 1)
polygon(c(myenvp$r[19:650], rev(myenvp$r[19:650])), 
        c((myenvp$lo/myenvp$mmean)[19:650], rev((myenvp$hi/myenvp$mmean)[19:650])),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = 1)
dev.off()


postscript(file=paste0(number, "_PCFStrSenResWholeSample.eps"))
plot(myenvp$r[19:650],(myenvp$obs/myenvp$mmean)[19:650],type = 'l',main = "PCF Sensitive and Resistant",
     xlab = expression(bold(paste("Stroma-Cell Distance"," ", "(", mu,"m",")"))), 
     cex.lab = 1.3, font.lab = 4,  ylim = c(0.0, 1.5), xlim = c(10, 320), col = "darkgreen", lwd = 3,
     ylab =""# expression(bold("g"["BrdU+,Str"])),
)
title(ylab=expression(bold("g"["Stroma,Sensitive/Resistant"])), line=2.5, cex.lab=1.4)
lines(myenv$r[19:650],(myenv$obs/myenv$mmean)[19:650],col = "red", lty = 1, lwd = 3)
# lines(x,(averObs)[15:PCF_r_length], col = rgb(0.3,0.2,1), lty = 1, lwd = 3)
# lines(x, y1, col = rgb(0.5,0.5,0.5), lty = 2, lwd = 1)
lines(myenvp$r[19:650],(myenvp$mmean/myenvp$mmean)[19:650], col = rgb(0.5,0.5,0.5), lty = 2, lwd = 1)
polygon(c(myenvp$r[19:650], rev(myenvp$r[19:650])), 
        c((myenvp$lo/myenvp$mmean)[19:650], rev((myenvp$hi/myenvp$mmean)[19:650])),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = 1)
dev.off()


XY <- superimpose(Pos=mypatternPos,Str=mypatternStr)
XYN <- superimpose(Neg=mypatternNeg,Str=mypatternStr)
PN <- superimpose(Pos=mypatternPos,Neg=mypatternNeg)
mysim <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Pos, Str = ppp2)
}

mysimN <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Neg, Str = ppp2)
}

mysimNP <- function(ppp1,ppp2){
  superimpose(PosR = split(rlabel(ppp1))$Neg, Pos = ppp2)
}

rJ <-seq(from = 0, to = 100, by = 0.5)
myenvJ <- envelope(XYN, fun = Jcross, correction=c("rs"),nsim=39, simulate=expression(mysimN(PN,mypatternStr)))
myenvpJ <- envelope(XY, fun = Jcross, correction=c("rs"),nsim=39, simulate=expression(mysim(PN,mypatternStr)))
myenvJPN <- envelope(PN, fun = Jcross, correction=c("rs"),nsim=39)#,r =rJ)#, simulate=expression(mysimNP(PN,mypatternPos)))

plot(myenvpJ)
plot(myenvJPN)

Z <- distmap(mypatternNeg)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,Z)
plot(myRho[1:200])
plot(myRhoP$Z, myRhoP$rho, xlim = c(0,300), ylim = c(0, 0.0008))
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
rLim <- 210


postscript(file=paste0(number, "_rhohatcombinedtwoAxesWSample.eps"))
par(mar = c(5, 4, 4, 4) + 0.3)              
plot(myRhoP$Z[1:rLim], myRhoP$rho[1:rLim],  main="Rho hat BrdU+ and BrdU-",
     xlab = expression(bold(paste("Distance from Stroma"," ", "(", mu,"m",")"))),
     #ylim= c(min(myRhoP$lo[1:rLim]),max(myRhoP$hi[1:rLim])),
     type='l', col= "darkgreen", lwd = 3,
     font.main=3, font.lab=2, cex.lab = 1.3,
     ylab="" #"BrdU- Density (# /unit area)",
)
title(ylab="Sensitive tumor Cells RSF",  line=2.5, cex.lab=1.3,font.lab=2, col.lab = "darkgreen")     
lines(myRhoP$Z[1:rLim], myRhoP$hi[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRhoP$Z[1:rLim], myRhoP$lo[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRhoP$Z[1:rLim], rev(myRhoP$Z[1:rLim])), c(myRhoP$lo[1:rLim], rev(myRhoP$hi[1:rLim]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRhoP$ave[20], lty = 2, lwd = 2, col = rgb(0,0.5,0))
par(new = TRUE)  
plot(myRho$Z[1:rLim], myRho$rho[1:rLim], type='l', col= "red", lwd = 3, 
     ylim= c(min(myRho$lo),max(myRho$hi)),   bty = "n",       
     axes = FALSE, xlab = "", ylab = "")
lines(myRho$Z[1:rLim], myRho$hi[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
lines(myRho$Z[1:rLim], myRho$lo[1:rLim], col = rgb(0.5,0.5,0.5), lty = 2)
polygon(c(myRho$Z[1:rLim], rev(myRho$Z[1:rLim])), c(myRho$lo[1:rLim], rev(myRho$hi[1:rLim]) ),
        col = rgb(0.5,0.5,0.5,0.3), lty = 0)
abline(h = myRho$ave[20], lty = 2, lwd = 2, col = rgb(0.9,0.2,0))
axis(side = 4)#, at = pretty(range(myRho$hi[1:rLim]))) 
mtext("Resistant Tumor Cells RSF",  side = 4, line=2.5, cex=1.3, font=2, col ="red" )
dev.off()

