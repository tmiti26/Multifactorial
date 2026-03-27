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
library(scales)

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
  
  xn <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Green')))]
  yn <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Green')))]
  xyn_inside<-inside.owin(xn,yn,all_owin)
  xxn <-xn[xyn_inside]
  yyn <- yn[xyn_inside]
  mypatternNeg <- ppp(xxn,yyn,window=all_owin)#
  numNeg <- numNeg + length(mypatternNeg$x)

  
  xs <- mydata$X[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  ys <- mydata$Y[apply(mydata, 1, function(r) any(r %in% c('Stroma')))]
  xys_inside<-inside.owin(xs,ys,all_owin)
  xxs <-xs[xys_inside]
  yys <- ys[xys_inside]
  mypatternStr <- ppp(xxs,yys,window=all_owin)#
  numStr <- numStr + length(mypatternStr$x)
  
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
  XPos_threshold_0_250 <- mypatternPos[nnn >= 0 & nnn <= 250]
  XPos_threshold_50_250 <- mypatternPos[nnn >= 50 & nnn <= 250]
  XPos_threshold_0_50 <- mypatternPos[nnn >= 0 & nnn <= 50]
  length(nnn)
  length(nnn[nnn>20])
  numPosStrk1 <- numPosStrk1 + length(nnn[nnn>20])
  mmm <- nncross(mypatternNeg, mypatternStrEdgeTot)$dist
  XNeg_threshold_0_250 <- mypatternNeg[mmm >= 0 & mmm <= 250]
  XNeg_threshold_50_250 <- mypatternNeg[mmm >= 50 & mmm <= 250]
  XNeg_threshold_0_50 <- mypatternNeg[mmm >= 0 & mmm <= 50]
  length(mmm)
  length(mmm[mmm>20])
  numNegStrk1 <- numNegStrk1 + length(mmm[mmm>20])
  
}


Z <- distmap(mypatternStr)
myDist <- as.im(Z)
myTest <- auc(mypatternPos, myDist, high = TRUE)
plot(myTestr <- roc(mypatternPos, myDist, high = TRUE))

b = quantile(myDist, prob = (0:5)/5) # a way to quantify stroma interdistance!!!!!!
bb = seq(from = 0, to = max(b), by = 25)
bVS <- c()
bVStotal <- c()
myPr <- (0:5)/5


# Build data frame from your vectors
df <- data.frame(
  Distance = b,
  CDF = myPr
)

# Plot
#cairo_ps("CummStrFract.eps", width = 6, height = 4, family = "Arial")

ggplot(df, aes(x = Distance, y = CDF)) +
  geom_point(shape = 21, fill = "white", color = "black", size = 2, stroke = 1) +
  labs(
    title = "Inter Stroma Distances",
    x = expression(bold("Distance")~"("~mu * "m"~")"),
    y = "Cummulative Fraction of Distances"
  ) +
  coord_cartesian(ylim = c(0, 1.0), xlim = c(0, 1500)) +
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


pdf(file= paste0(i,"_Plot.pdf"))
plot(mypatternNeg, type = "p", cex = 0.00001, col = "red")
plot(mypatternPos, type = "p", cex = 0.02, col = "green", add = TRUE)
plot(mypatternStr, type = "p", cex = 0.02, col = "blue", add = TRUE)
dev.off()

# Convert point patterns to data frames
df_neg <- as.data.frame(mypatternNeg)
df_pos <- as.data.frame(mypatternPos)
df_str <- as.data.frame(mypatternStr)
# Add a group column to identify each pattern
df_neg$group <- "Resistant"
df_pos$group <- "Sensitive"
df_str$group <- "Stromal"
# Combine all into one data frame
df_all <- rbind(df_neg, df_pos, df_str)
# Assign sizes and colors
df_all$size <- ifelse(df_all$group == "Resistant", 0.00001, 0.00002)
df_all$color <- ifelse(df_all$group == "Resistanr", "red",
                       ifelse(df_all$group == "Sensitive", "green", "blue"))
# Plot
ggplot(df_all, aes(x = x, y = y)) +
  geom_point(aes(color = group, size = size)) +
  scale_color_manual(values = c("Resistant" = "red", "Sensitive" = "green", "Stromal" = "blue")) +
  scale_size_identity() +
  coord_fixed() +
  theme_minimal()
ggsave(filename = "334_SamplePlot.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

nn <- nncross(mypatternPos, mypatternStr)$dist
mm <- nncross(mypatternNeg, mypatternStr)$dist
for(jj in 1:length(nn)){
  nncrosss <- c(nncrosss, nn[jj])}
for(hh in 1:length(mm)){
  nncrossr <- c(nncrossr, mm[hh])}

xRho <- 202
nn_0_50SR <- nncross(XPos_threshold_0_50, XNeg_threshold_0_50)$dist
plot(density(nn_0_50SR))
nn_0_50SStr <- nncross(XPos_threshold_0_50, mypatternStr)$dist
plot(density(nn_0_50SStr))
nn_0_50RStr <- nncross(XNeg_threshold_0_50, mypatternStr)$dist
plot(density(nn_0_50RStr))


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


#cairo_ps("334_Dist0_50.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("RS" = "black", "SStr" = "darkgreen", "RStr" = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, xRho), ylim = c(0, 0.035)) +
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

ggsave(filename = "334_Dist0_50.eps", 
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


#cairo_ps("334_Dist50_250.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("RS" = "black", "SStr" = "darkgreen", "RStr" = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, xRho), ylim = c(0, 0.035)) +
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

ggsave(filename = "334_Dist50_250.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
################


nn_0_250SR <- nncross(XPos_threshold_0_250, XNeg_threshold_0_250)$dist
plot(density(nn_0_250SR))
nn_0_250SStr <- nncross(XPos_threshold_0_250, mypatternStr)$dist
plot(density(nn_0_250SStr))
nn_0_250RStr <- nncross(XNeg_threshold_0_250, mypatternStr)$dist
plot(density(nn_0_250RStr))

#######ggplot
# Compute density objects
dens_RS   <- density(nn_0_250SR)
dens_SStr <- density(nn_0_250SStr)
dens_RStr <- density(nn_0_250RStr)

# Combine into one long-format data frame
df <- data.frame(
  x = c(dens_RS$x, dens_SStr$x, dens_RStr$x),
  y = c(dens_RS$y, dens_SStr$y, dens_RStr$y),
  group = rep(c("Resistant - Sensitive", "Sensitive - Str.", "Resistant - Str."), each = length(dens_RS$x))
)


#cairo_ps("334_Dist50_500.eps", width = 8, height = 6, family = "Arial")

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Resistant - Sensitive" = "black", "Sensitive - Str." = "darkgreen",
                                "Resistant - Str." = "red")) +
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density",
    color = ""  # Legend title
  ) +
  coord_cartesian(xlim = c(0, xRho), ylim = c(0, 0.035)) +
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

ggsave(filename = "334_Dist0_250.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
################




nn <- density(nncrosss)
mm <- density(nncrossr)
xmaxs <- which.max(density(nncrosss)$y)
density(nncrosss)$x[xmaxs]
xmaxr <- which.max(density(nncrossr)$y)
density(nncrossr)$x[xmaxr]

pdf(file=paste0(number,"_DistToStromaWholeSample.pdf"),width = 8, height = 8)
plot(density(nncrosss), col = "darkgreen", lwd = 3, xlim = c(-0, xRho),
     xlab = expression(bold(paste("r"," ", "(", mu,"m",")"))),
     ylab= "Sensitive/Resistant Density", font.lab=2)
lines(density(nncrossr), col = "red", lwd = 3)
polygon(nn, density = -1, col = rgb(0.4,0.6,0.2, 0.5))  
polygon(mm, density = -1, col = rgb(0.9, 0, 0.1, 0.4))
abline(v = density(nncrosss)$x[xmaxs], col = "darkgreen", lty = 1, lwd = 2)
abline(v = density(nncrossr)$x[xmaxr], col = "red", lty = 1, lwd = 2)
dev.off()



################### ggplot
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

# Plot
#cairo_ps("334_DistToStromaWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot() +
  
  # Density lines
  geom_line(data = df_density, aes(x = x, y = y, color = group), size = 1.3) +
  scale_color_manual(values = c("Sensitive" = "darkgreen", "Resistant" = "red")) +
  
  labs(
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Density"
  ) +
  coord_cartesian(xlim = c(0, xRho)) +
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

f <- ecdf(nncrosss)
ff <- ecdf(nncrossr)
myKS <- ks.test(nncrosss, ff)

##########ggplot
# Convert step function to data frame
df_f <- data.frame(x = knots(f), y = f(knots(f)), group = "Sensitive")
df_ff <- data.frame(x = knots(ff), y = ff(knots(ff)), group = "Resistant")
df_all <- rbind(df_f, df_ff)

#cairo_ps("334_CumKSWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot(df_all, aes(x = x, y = y, color = group)) +
  geom_step(size = 1.3) +  # Use step to mimic base R's ecdf style
  scale_color_manual(values = c("Sensitive" = "darkgreen", "Resistant" = "red")) +
  labs(
    title = "NN Cumulative Density Function",
    x = expression(bold(paste("r", " ", "(", mu, "m", ")"))),
    y = "Sensitive/Resistant Cumulative Density"
  ) +
  coord_cartesian(xlim = c(0, xRho)) +
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


################ggplot
# Get min and max for y-values
Z <- distmap(mypatternStr)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,Z)
plot(myRho[1:100])
plot(myRhoP[1:100])
rLim <- 210
rRho <-295
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

#cairo_ps("334_rhohatCombinedtwoAxesWholeSample.eps", width = 8, height = 6, family = "Arial")
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
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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
rLim <- 210
rRho <-300
range_sensitive <- range(myRhoP$rho[1:rRho])
range_resistant <- range(myRho$rho[1:rRho])

df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
#cairo_ps("334_rhohatSRWholeSample0_50.eps", width = 8, height = 6, family = "Arial")
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
  scale_y_continuous(labels = scientific) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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
rLim <- 210
rRho <-295
range_sensitive <- range(myRhoP$rho[1:rRho])
range_resistant <- range(myRho$rho[1:rRho])
df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
#cairo_ps("334_rhohatSRWholeSample50_250.eps", width = 8, height = 6, family = "Arial")
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
  scale_y_continuous(labels = scientific) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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
myRhoP <- rhohat(XPos_threshold_0_250,Z)
rLim <- 210
rRho <-295
range_sensitive <- range(myRhoP$rho[1:rRho])
range_resistant <- range(myRho$rho[1:rRho])
df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
#cairo_ps("334_rhohatSRWholeSample0_250.eps", width = 8, height = 6, family = "Arial")
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
  scale_y_continuous(labels = scientific) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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

ggsave(filename = "334_rhohatSRWholeSample0_250.eps",
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

#################

#######ggplot
Z <- distmap(mypatternNeg)
myRhoP <- rhohat(mypatternPos,Z)
rLim <- 210
rRho <-295
range_sensitive <- range(myRhoP$rho[1:rRho])
range_resistant <- range(myRho$rho[1:rRho])
df <- data.frame(
  Z = myRhoP$Z,
  rho = myRhoP$rho,
  lo = myRhoP$lo,
  hi = myRhoP$hi
)
#cairo_ps("334_rhohatSRWholeSample0_250.eps", width = 8, height = 6, family = "Arial")
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
  scale_y_continuous(labels = scientific) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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

ggsave(filename = "334_rhohatSRWholeSample.eps",
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")

#################
#######ggplot combined positive troma and positive resistant

# Get min and max for y-values
Z <- distmap(mypatternStr)
ZZ <- distmap(mypatternNeg)
myRho <- rhohat(mypatternPos,Z)
myRhoP <- rhohat(mypatternPos,ZZ)
plot(myRho[1:100])
plot(myRhoP[1:100])
rLim <- 210
rRho <-295
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
  group = rep(c("Stroma", "Resistant"), each = rRho)
)

ref_lines <- data.frame(
  yintercept = c(myRhoP$ave[20], scaled_ave),
  group = c("Stroma", "Resistant")
)

#cairo_ps("334_rhohatCombinedtwoAxesWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot(df, aes(x = Z, y = rho, group = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  geom_line(aes(color = group), size = 1.2) +
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  geom_hline(data = ref_lines, aes(yintercept = yintercept, color = group),
             linetype = "dashed", size = 1) +
  
  scale_color_manual(values = c("Stroma" = "darkgreen", "Resistant" = "black")) +
  scale_y_continuous(
    name = "RSF Sensiive/Stroma",
    labels = scales::scientific,
    sec.axis = sec_axis(
      ~ (. - offset) / scale_factor,
      name = "RSF Sensitive/Resistant Cells",
      labels = scales::scientific
    )
  ) +
  
  labs(
    title = "Rho hat Stroma and Resistant Tumor Cells",
    x = expression(bold(paste("Distance (", mu, "m)")))
  ) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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
    axis.title.y.right = element_text(face = "bold", size = 13, color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
ggsave(filename = "334_rhohatCombinedSStrStwoAxesWholeSample.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")


#######ggplot combined resistant stroma and positive resistant

# Get min and max for y-values
Z <- distmap(mypatternStr)
ZZ <- distmap(mypatternNeg)
myRho <- rhohat(mypatternNeg,Z)
myRhoP <- rhohat(mypatternPos,ZZ)
plot(myRho[1:100])
plot(myRhoP[1:100])
rLim <- 210
rRho <-295
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
  group = rep(c("Stroma", "Resistant"), each = rRho)
)

ref_lines <- data.frame(
  yintercept = c(myRhoP$ave[20], scaled_ave),
  group = c("Stroma", "Resistant")
)

#cairo_ps("334_rhohatCombinedtwoAxesWholeSample.eps", width = 8, height = 6, family = "Arial")
ggplot(df, aes(x = Z, y = rho, group = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey50", alpha = 0.3) +
  geom_line(aes(color = group), size = 1.2) +
  geom_line(aes(y = hi), color = "grey50", linetype = "dashed") +
  geom_line(aes(y = lo), color = "grey50", linetype = "dashed") +
  geom_hline(data = ref_lines, aes(yintercept = yintercept, color = group),
             linetype = "dashed", size = 1) +
  
  scale_color_manual(values = c("Stroma" = "red", "Resistant" = "black")) +
  scale_y_continuous(
    name = "RSF Resistant/Stroma",
    labels = scales::scientific,
    sec.axis = sec_axis(
      ~ (. - offset) / scale_factor,
      name = "RSF Sensitive/Resistant Cells",
      labels = scales::scientific
    )
  ) +
  
  labs(
    title = "Rho hat Stroma and Resistant Tumor Cells",
    x = expression(bold(paste("Distance (", mu, "m)")))
  ) +
  coord_cartesian(xlim = c(0, xRho)) +
  scale_x_continuous(breaks = seq(0, max(df_all$x), by = 50)) +  
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
    axis.title.y.left = element_text(face = "bold", size = 13, color = "red"),
    axis.title.y.right = element_text(face = "bold", size = 13, color = "black"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
ggsave(filename = "334_rhohatCombinedRStrRtwoAxesWholeSample.eps", 
       device = cairo_ps, 
       dpi = 1200, 
       width = 8,
       height = 6, 
       family = "Arial")
