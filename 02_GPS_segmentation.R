#This code downloads honey buzzard GPS data and classifies them into flight types
#follows Martina's code: https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS1_soaring_skills/script1_GPSsegmentation_goldenEagles_newMarch2023.R
#Elham Nourani PhD.
#Mar 16. 2023. Konstanz, DE.a

library(move)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(plyr)
library(data.table)
library(doParallel)
detectCores()
doParallel::registerDoParallel(detectCores()-2) 
#dir.create("prepped")
#library(data.table) #for Martina's rbindlist code
#download gps and IMU from movebank. focus on one individual: D324-512
#setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")
setwd("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#check for sensors available. Acc = 2365683; Mag = 77740402; Orientation: 819073350
getMovebankSensors(EHB_FN_id,login = creds)

#mv_ls <- getMovebankData(EHB_FN_id,  animalName = c("D324_512", "D320_475"), removeDuplicatedTimestamps = T, login = creds)
#download data for all individuals
mv_ls <- getMovebankData(EHB_FN_id, removeDuplicatedTimestamps = T, login = creds)

saveRDS(mv_ls,"/home/enourani/Desktop/HB_data_Apr23/all_GPS_Apr4.rds")

mv_ls <- split(mv_ls)

#keep only move objects (individuals) with more than 30 locations (needed for segmentation)
mv_ls <- mv_ls[sapply(mv_ls, n.locs) >= 30]

dir.create("GPS_seg_Apr23/classifiedData/")

(b <- Sys.time())

lapply(mv_ls, function(mv){
  
  # STEP 1: calculate track variables -------------------------------------------------------------------------
  animalID <- mv@idData$local_identifier
  mv$timelag.sec <- c(NA, timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA, (mv$height_above_ellipsoid[-1] - mv$height_above_ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA, distance(mv))
  mv$gr.speed <- c(NA, speed(mv))
  
  #STEP 2: # Segmentation on bursts of continuous 1 sec resolution data ---------------------------------------
  minResol <- 2 # 1 to max 2 sec timelag
  minBurstDuration <- 30 # we want bursts of at least 30 secs
  
  swV <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
  swT <- 12 #smoothing window of 29 seconds for thermalling behaviour (according to Rolf/Bart's paper) (we changed it to smooth window of 21 sec)
  circlDegrees <- 250 #degrees of rotation to be achieved in the time defined by swT*2
  minBehavDuration <- 5 #minimum duration in seconds of a specific behaviour, when less than this and if in between two segments of a different behaviour it will be incorporated in the previous and follwoing segment 
  minThermalDuration <- 30 #minimum duration for a circling event to be considered as thermalling
  
  # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
  mv$burstID <- c(0, cumsum(mv$timelag.sec[2:n.locs(mv)] > minResol))  #from row 2 (since the first is NA)
  # with table we can count all locations with the same ID (each one is a separate burst) and keep those high resolution bursts that have at least a certain number of locations (= minBurstDuration)
  burstDuration <- as.data.frame(table(mv$burstID))
  burstsToKeep <- burstDuration$Var1[which(burstDuration$Freq >= minBurstDuration)]
  # use those to subset the move obj and keep only the high resolution bursts
  HRmv <- mv[which(mv$burstID %in% burstsToKeep),]
  if(nrow(HRmv)>0){
    HRdf_bursts <- as.data.frame(HRmv)
    # Remove unnecessary columns
    HRdf_bursts <- HRdf_bursts[,-grep("mag|orientation|coords|timestamps|start.timestamp|optional|import|visible|algorithm|battery|decoding|accuracy|manually|activity|checksum|acceleration",
                                      colnames(HRdf_bursts))]
    # Split each individual dataframe by burst ID
    burst_ls_corr <- split(HRdf_bursts, HRdf_bursts$burstID)
    # Keep only bursts with minBurstDuration (30 of smoothing window will be NA) 
    burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= minBurstDuration)]
    # Compute smoothed turning angle separately for each burst
    HRdf <- llply(burst_ls_corr_sub, function(b){ # alternatively run in parallel
        b$vertSpeed_smooth <- NA
        b$turnAngle_smooth <- NA
        for(i in (swV+1):(nrow(b)-swV)){
          b$vertSpeed_smooth[i] <- mean(b$vert.speed[(i-swV):(i+swV)], na.rm=T)}
        for(i in (swT+1):(nrow(b)-swT)){
          b$turnAngle_smooth[i] <- max(abs(cumsum(b$turn.angle[(i-swT):(i+swT)])))}
        return(b) # return df with smoothed variables
      }, .parallel = T) %>% 
        reduce(rbind) %>% 
        as.data.frame()
    
    # Classify soaring only based on vertical speed
    HRdf <- HRdf[complete.cases(HRdf$vertSpeed_smooth),]
    kmeanV <- kmeans(HRdf$vertSpeed_smooth, 2)   #Get the two clusters
    soarId <- which.max(aggregate(HRdf$vertSpeed_smooth~kmeanV$cluster, FUN=mean)[,2]) # which one is the soaring one?
    soarClust <- rep("glide", length(kmeanV$cluster))
    soarClust[which(kmeanV$cluster==soarId)] <- "soar"
    HRdf$soarClust <- factor(soarClust, levels=c("soar","glide"))  
    # Now classify thermalling only based on turning angle (cumulated to a 25 s time window in previous step)
    HRdf$flightClust <- "other"
    HRdf$flightClust[which(HRdf$soarClust=="soar" & HRdf$turnAngle_smooth >= circlDegrees)] <- "circular soaring" #complete 150 degrees in 15 sec
    HRdf$flightClust[which(HRdf$soarClust=="soar" & HRdf$flightClust != "circular soaring")] <- "linear soaring"
    HRdf$flightClust[which(HRdf$soarClust=="glide")] <- "gliding"
    HRdf$flightClust <- factor(HRdf$flightClust, levels=c("circular soaring","linear soaring","gliding","other"))
    
    # Add some steps of smoothing based on duration of behaviours:
    burst_ls_class <- split(HRdf, HRdf$burstID)
    burst_ls_class_smooth <- llply(burst_ls_class, function(b){
      # We assign a unique ID to each consecutive flight segment based on a rule (the ID increases if the class of each obs is different from the previous one)
      print(unique(b$burstID))
      b <- b[order(b$timestamp),]
      b$flightNum <- c(0, cumsum(b$flightClust[-1] != b$flightClust[-nrow(b)]))
      
      # we calculate the duration and unique class of each behavioural segment
      behavDuration <- merge(aggregate(timelag.sec~flightNum, data=b, FUN=sum), 
                             aggregate(flightClust~flightNum, data=b, FUN=unique), by="flightNum")
      # we create a new category ID, which is the same as the original one, 
      # unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
      behavDuration$flightNum_smooth <- behavDuration$flightNum
      behavDuration$flightClust_smooth <- behavDuration$flightClust
      if(nrow(behavDuration)>2){
        for(i in 2:(nrow(behavDuration)-1)){
          if(behavDuration$timelag.sec[i] <= minBehavDuration & behavDuration$flightClust_smooth[i-1] == behavDuration$flightClust_smooth[i+1]){
            behavDuration$flightNum_smooth[c(i, i+1)] <- behavDuration$flightNum_smooth[i-1]
            behavDuration$flightClust_smooth[c(i, i+1)] <- behavDuration$flightClust_smooth[i-1]
          }}
      }
      if(nrow(behavDuration)==2){
        if(any(behavDuration$timelag.sec <= minBehavDuration)){
          longestBehav <- which.max(behavDuration$timelag.sec)
          behavDuration$flightNum_smooth <- behavDuration$flightNum_smooth[longestBehav]
          behavDuration$flightClust_smooth <- behavDuration$flightClust_smooth[longestBehav]
        }
      }
      b <- merge(b, behavDuration[c("flightNum","flightNum_smooth","flightClust_smooth")], by="flightNum", all.x=T)
      
      # recalculate segment duration based on smoothed classification and reclassify as linear soaring all circling that lasts < 30 seconds
      behavDuration_smooth <- merge(aggregate(timelag.sec~flightNum_smooth, data=b, FUN=sum), 
                                    aggregate(flightClust_smooth~flightNum_smooth, data=b, FUN=unique), by="flightNum_smooth")
      behavDuration_smooth$flightNum_smooth2 <- behavDuration_smooth$flightNum_smooth
      behavDuration_smooth$flightClust_smooth2 <- behavDuration_smooth$flightClust_smooth
      if(nrow(behavDuration_smooth)>2){
        for(i in 2:(nrow(behavDuration_smooth)-1)){
          if(behavDuration_smooth$flightClust_smooth2[i] == "circular soaring"){
            if(behavDuration_smooth$timelag.sec[i] <= minThermalDuration & behavDuration_smooth$flightClust_smooth2[i-1] == "linear soaring" & behavDuration_smooth$flightClust_smooth2[i+1] == "linear soaring"){
              behavDuration_smooth$flightNum_smooth2[c(i, i+1)] <- behavDuration_smooth$flightNum_smooth2[i-1]
              behavDuration_smooth$flightClust_smooth2[c(i, i+1)] <- behavDuration_smooth$flightClust_smooth2[i-1]
            }}}
      }
      b <- merge(b, behavDuration_smooth[c("flightNum_smooth","flightNum_smooth2","flightClust_smooth2")], by="flightNum_smooth", all.x=T)
      # finally check the classification of the time window at the start and end of the track
      # these first and last points can only be classified as either gliding or linear, as their classification was only based on vertical speed but not turning angle
      # so if at the start or end there is linear soaring, but they are preceded or followed by circluar soaring, they become circular
      behavDuration_smooth2 <- merge(aggregate(timelag.sec~flightNum_smooth2, data=b, FUN=sum), 
                                     aggregate(flightClust_smooth2~flightNum_smooth2, data=b, FUN=unique), by="flightNum_smooth2")
      behavDuration_smooth2$flightNum_smooth3 <- behavDuration_smooth2$flightNum_smooth2
      behavDuration_smooth2$flightClust_smooth3 <- behavDuration_smooth2$flightClust_smooth2
      if(nrow(behavDuration_smooth2) >= 2){
        if(behavDuration_smooth2$timelag.sec[1] <= swT & 
           behavDuration_smooth2$flightClust_smooth2[1] == "linear soaring" & behavDuration_smooth2$flightClust_smooth2[2] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[1] <- behavDuration_smooth2$flightNum_smooth3[2]
          behavDuration_smooth2$flightClust_smooth3[1] <- "circular soaring"
        }
        if(behavDuration_smooth2$timelag.sec[nrow(behavDuration_smooth2)] <= swT & 
           behavDuration_smooth2$flightClust_smooth2[nrow(behavDuration_smooth2)] == "linear soaring" & behavDuration_smooth2$flightClust_smooth2[nrow(behavDuration_smooth2)-1] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)] <- behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)-1]
          behavDuration_smooth2$flightClust_smooth3[nrow(behavDuration_smooth2)] <- "circular soaring"
        }}
      # merge with burst
      b <- merge(b, behavDuration_smooth2[c("flightNum_smooth2","flightNum_smooth3","flightClust_smooth3")], by="flightNum_smooth2", all.x=T)
      
      # Assign unique ID to the behavioural segment based on the final smoothest classification
      b$track_flight_id <- paste0(unique(b$individual_id),"_",unique(b$burstID),"_segm_",b$flightNum_smooth3) 
      
      return(b) #return each classified and smoothed burst to a list
    }, .parallel=T) # to make it run in parallel use llply (instead of lapply) with .parallel=T
    
    # Rbind all bursts and save classified and smoothed dataframe per individual
    HRdf_smooth <- as.data.frame(rbindlist(burst_ls_class_smooth))
    save(HRdf_smooth, file = paste0("GPS_seg_Apr23/classifiedData/animal_",animalID,"_classifiedBursts_df.rdata"))
  }
})

Sys.time() - b #1.7 hrs

#STEP 3: plotting ------------------------------------------------------------------------
  
  burst_ls_class_smooth <- split(HRdf_smooth, HRdf_smooth$burstID)
  randomBursts <- sample(1:length(burst_ls_class_smooth), 100)
  
  lapply(burst_ls_class_smooth[randomBursts], function(b){
    #b=burst_ls_class_smooth[["7006"]]
    
  print(unique(b$burstID))
  animalID <- unique(b$local_identifier)
    
  b <- b[order(b$timestamp),]
  # cbind(as.character(b$flightClust),as.character(b$flightClust_smooth),as.character(b$flightClust_smooth2), as.character(b$flightClust_smooth3))
  
  # calculate aspect ratio for plot
  rangeLong <- max(b$location_long)-min(b$location_long)
  rangeLat <- max(b$location_lat)-min(b$location_lat)
  
  # Plot results
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),".png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
  dev.off()
  
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth.png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
  dev.off()
  
  # this 3D plots are only useful when interactive, exporting them doesn not make much sense
  # plot3d(b[,c("location_long","location_lat","height_above_ellipsoid")], type="l", col="darkgrey")
  # points3d(b[,c("location_long","location_lat","height_above_ellipsoid")], col=c("red","darkgreen","blue","grey")[b$flightClust_smooth3], size=5)
  # aspect3d(x=rangeLong/rangeLat, y=1, z=1)
  # snapshot3d(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth3D.png"), width = 600, height = 600)
  })  
  
