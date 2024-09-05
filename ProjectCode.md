Combining graph theory and spatially-explicit, individual-based models to improve
invasive species control strategies at a regional scale
================

## Table of Contents

-   Introduction
-   Code for generating simulated landscapes
-   Code for combining census files
-   Code for analysis of derived HexSim output files

## Introduction and Caveats

This markdown includes code used to collate, shape, and analyze HexSim
output data associated with the manuscript: Combining graph theory and spatially-explicit,
individual-based models to improve invasive species control strategies at a regional scale
This manuscript focuses on exploring optimal management strategies of the American bullfrog, *Lithobates catesbaeinus*,
informed by spatial data from the Huachuca Mountains and Canelo Hills in
southeastern Arizona.

The first section details how we simulated multiple landscapes for use in the IBM.
As a result of the need of multiple replicates across landscape
iterations and for each sensitivity analysis and for each simulated
management scenario, there was a bounty of data needed to be shaped and
combined to allow for downstream analyses in the final sections. In the second section, we
present the code to replicate the analyses reported in the manuscript.


While we present these code for transparency of methodology and
potential replication, newer version of HexSim have added increased
functionality rendering some of the workarounds (which were required the
carpentry of the data we analyzed derived from the individually based
model and its iteration) not necessarily needed. Read the annotations to
better understand different aspects of the code and what was performed.
We would recommend these scripts being embedded into an R project for
convenience of directories working to your benefit.

Preferred Citation:

Drake, J.C., O'Malley, G., Kraft, J., Mims, M.C. 2024. Combining graph theory and spatially-explicit,
individual-based models to improve invasive species control strategies at a regional scale. *Landscape Ecology* (In Revision).

##

```r
# Drake et al. 2024: sIBMS for LICA control 
# Appendix S3 Pt1
# Simulated landscape generation 
# For SMARTSIM associated manuscript produced for Dr. Meryl Mims Lab
# March 2024
# Copyright J.Drake 2024 - noncommercial nonprofit use permitted 

library("terra")
library("dplyr")
library("tidyr")
library("sf")
library("ggplot2")
library("tidyterra")



set.seed(42) # set seed for reproducibility 

#### Create landscape elements ----
  
  #* Create uniform landscape raster ----
    # Rectangle extent = 30km^2 
    # Cell size = 30 meter
    # Cell value = 1 ( or 0? )

    # for 30mX30m grid nrow&ncol=1000
    # for 5mX5m grid = 6000 rows
    # fpr 1x1m grid = 30,000 rows cols
    

n.col <- 1000 # if needed can increase the data to 6000 to make landscape cell size smaller
n.row <- 1000 # if needed can increase the data to 6000 to make landscape cell size smaller
x.min <- 0
y.min <- 0
x.max <- 30000
y.max <- 30000
vals <- 1
background <- -999 # for mode to import waters correctly, we will use
                   # a special background for the water patch generation
                   # specifically, 999, this allows you to identify in hexsim 
                   # this cell designation to define as unscored matrix

r <- rast(ncol=n.col,
          nrow=n.row,
          xmin=x.min,
          xmax=x.max,
          ymin=y.min,
          ymax=y.max,
          #crs=,
          #extent=,
          #resolution=,
          vals=vals,
            )


  #* Create management area in center ----
# define the center of the map
center.pt <- matrix(data=NA, nrow=1, ncol=2) 
center.pt[,1] <- x.max/2
center.pt[,2] <- y.max/2
center.pt <- vect(center.pt)
  
plot(r, asp=1)
plot(center.pt, add=TRUE)

# radius = sqrt(A/pi)
radius <- sqrt((0.15*x.max*y.max)/pi) # radius set to encapsulate 15% of the total landscape area

# create the buffer
center.bf <- buffer(center.pt, radius)
MgmtArea <- mask(x=r,
                 mask=center.bf
                 )

MgmtArea
plot(r)
plot(MgmtArea, add=TRUE, col="red")
plot(center.pt, add=TRUE, col="black")


  #* Create static waters ----
    # Uniform in size = 30m (1 cell)
    # Number of waters = Approximating study area ~600 ponds
      # 600*0.45 = 270 waters
      # Of those ponds with hydroperiod info available (n=155)
      # 62.58065% of the those had hydro index is greater than 0.6
      # 0.63*270 = 170.1
    
# Distribution of waters = uniform random Northings & Eastings

num.patches <- 170 # total number of patches
min.border <- 100 # to prevent waters from falling outside boundary when imported to hexsim
max.border <- 29900 # to prevent waters from falling outside boundary when imported to hexsim
n.iter <- 100 # number of landscapes wanted, we generated 100, but ended up using only the first  10.

    # create bucket for generated data

simpatch <- array(NA, dim=c(num.patches, 2, n.iter)) # 100 patches (rows), 2 fields xy, and 100 iterations

    # for loop for generating landscape 

for(i in 1:n.iter) {
  simpatch[,1,i] <- round(runif(num.patches, min=min.border, max=max.border), 0)
  simpatch[,2,i] <- round(runif(num.patches, min=min.border, max=max.border), 0)
}
simpatch

if(1==2){ # turn 'on' for generation of files for hexsim import
  

  for(i in 1:dim(simpatch)[3]){
    # this is used to generate files for import to hexsim to show locaiton of all habitat patches
    temp.patches <- 
      rasterize(x=simpatch[,,i],
                y=r,
                values=1,
                background=-999) # special background value for use with mode import in HexSIM
    writeRaster(temp.patches,
                filename = paste0("SimulationLandscape/WaterSites/AllOne/Patches",i,".asc" ),
                NAflag = -999,
                overwrite = TRUE 
                )
    # this is used for import to generate individual patch id's
    temp.ind.patches <- 
      rasterize(x=simpatch[,,i],
                y=r,
                values=1:num.patches,
                background=-999) # special background value for use with mode import in HexSIM
    writeRaster(temp.ind.patches,
                filename = paste0("SimulationLandscape/WaterSites/IndNumbers/IndPatches",i,".asc" ),
                NAflag = -999,
                overwrite = TRUE 
    )

    
  }

}



# Management Scenarios ----------------------------------
n.col <- 1000
n.row <- 1000
x.min <- 0
y.min <- 0
x.max <- 30000
y.max <- 30000
vals <- 1
background <- -999 # for 'mode' variant of import for our model, we will use
# a special background for the water patch generation
# specifically, 999

r <- rast(ncol=n.col,
          nrow=n.row,
          xmin=x.min,
          xmax=x.max,
          ymin=y.min,
          ymax=y.max,
          #crs=,
          #extent=,
          #resolution=,
          vals=vals,
)

  # We will use a 5%, 20%, and 50% management prioritization schema 
  # where we are focusing effort on the top 5%, top 20%, or top 50% of waters
    # 5% = 170*0.05 = 8.5 = 9
    # 20% = 170*0.2 = 34
    # 50% = 170*0.5 = 85


  #* Connectivity prioritized management treatments----
    # constructing the graph and deciding edge rules for node connections

      # for now edge rules is 5 kilometer connection, and I am going to think about it

    # will use eigenvalue centrality 
    # will use betweenness centrality
library(igraph)
  #* Connectivity metrics for loop ----

dim(simpatch[,,])

  # the set up info

num.patches <- 170
  # choose which control effort level we will run 
#perc.ctrl.patches <- 0.05 
#perc.ctrl.patches <- 0.2   
perc.ctrl.patches <- 0.5 

n <- round(num.patches*perc.ctrl.patches, 0)
n

  # the for loop

for (i in 1:dim(simpatch)[3]){

temp.patches <- simpatch[,,i]
d.pts <- dist(temp.patches)
dm.pts <- as.matrix(d.pts)
# set max edge distance using hard cutoff
dm.pts[dm.pts > 5000] <- 0
g.pts <- graph.adjacency(dm.pts,
                         mode = c("undirected"),
                         weighted = TRUE)
 
#plot(g.pts)
#associate with point locations
temp.patches <- as.data.frame(temp.patches)
temp.patches$btwnCentrality <- as.matrix(betweenness(g.pts))
temp.patches$eigenCentrality <- as.matrix(eigen_centrality(g.pts)$vector)
temp.patches

sorted.eig.rank=order(eigen_centrality(g.pts)$vector, decreasing = TRUE)
x<- temp.patches[,1][sorted.eig.rank[1:n]]
y<- temp.patches[,2][sorted.eig.rank[1:n]]
#temp.eig.patches <- buffer(vect(cbind(x,y)), 30)
temp.eig.patches <- 
  rasterize(x=buffer(vect(cbind(x,y)), 30),
          y=r,
          values=1,
          background=background)
writeRaster(temp.eig.patches, # choose the file path for correct control effort
            #filename = paste0("SimulationLandscape/ControlAreas/Eigenvalue/eigTop5perc",i,".asc" ),
            #filename = paste0("SimulationLandscape/ControlAreas/Eigenvalue/eigTop20perc",i,".asc" ),
            filename = paste0("SimulationLandscape/ControlAreas/Eigenvalue/eigTop50perc",i,".asc" ),
            NAflag = -999,
            overwrite = TRUE 
            )

sorted.btw.rank <- order(betweenness(g.pts), decreasing = TRUE) 
x <- temp.patches[,1][sorted.btw.rank[1:n]]
y <- temp.patches[,2][sorted.btw.rank[1:n]]
#temp.btwn.patches <- cbind(x,y)
temp.btwn.patches <- 
  rasterize(x=buffer(vect(cbind(x,y)), 30),
            y=r,
            values=1,
            background=background)
writeRaster(temp.btwn.patches, # choose the file path for correct control effort
            #filename = paste0("SimulationLandscape/ControlAreas/Betweenness/btwnTop5perc",i,".asc" ),
            #filename = paste0("SimulationLandscape/ControlAreas/Betweenness/btwnTop20perc",i,".asc" ),
            filename = paste0("SimulationLandscape/ControlAreas/Betweenness/btwnTop50perc",i,".asc" ),
            NAflag = -999,
            overwrite = TRUE 
)



}
#}


  #* Distance based measures for management treatments ----

library(spatstat)
library(matrixStats) # for rowMins function

  # the set up info
dim(simpatch[,,])

num.patches <- 170
  
  # choose which control level we will run

#perc.ctrl.patches <- 0.05 
#perc.ctrl.patches <- 0.2   
perc.ctrl.patches <- 0.5 

n <- round(num.patches*perc.ctrl.patches, 0)
n

for (i in 1:dim(simpatch)[3]){
  d.pts <- dist(simpatch[,,i])
  d.m.pts <- as.matrix(d.pts)
  diag(d.m.pts) <-  NA
  closestPatches <- simpatch[,,i][order(rowMins(d.m.pts, na.rm=TRUE, value=TRUE), decreasing = FALSE)[1:n],] 
  temp.closest.patches <- 
    rasterize(x=buffer(vect(closestPatches), 30),
              y=r,
              values=1,
              background=background)
  writeRaster(temp.closest.patches, # choose the file path for correct control effort
              #filename = paste0("SimulationLandscape/ControlAreas/ClosestPatches/closestTop5perc",i,".asc" ),
              #filename = paste0("SimulationLandscape/ControlAreas/ClosestPatches/closestTop20perc",i,".asc" ),
              filename = paste0("SimulationLandscape/ControlAreas/ClosestPatches/closestTop50perc",i,".asc" ),
              NAflag = -999,
              overwrite = TRUE 
  )
  
}

  #* Buffer distance around mgmt area ----



  # set number of waters to find
# 5% = 170*0.05 = 8.5 = 8 in R
# 20% = 170*0.2 = 34
# 50% = 170*0.5 = 85

perc.ctrl.patches <- 0.05 
#perc.ctrl.patches <- 0.2   
#perc.ctrl.patches <- 0.5 
n <- round(num.patches*perc.ctrl.patches, 0)
n
  # buffer around management area by 100 meter increment

# below is a little function to automate a search window around our protected are
# to generate the bespoke buffer based on the amount of patches defined by our 
# level of control effort

center.bf <- st_buffer(st_as_sf(center.pt), radius)

buffered.ctrl.patches <- function(priority.area=center.bf,
                                  waters=st_multipoint(pts),
                                  total.patches=170,
                                  percent.ctrl=0.05,
                                  increment=10,
                                  min.ring=10,
                                  max.ring=10000)
{
  #priority.area = area around you construct buffer
  #total.patches = total number of patches in landcape
  #percent.ctrl = percent of patches you want controlled
  #increment = increment to increase search window in landscape units (meters)
  
  # necessary libraries
  require(sf)
  require(terra)
  
  # Create container for results (control patches found in buffer ring)
  control.patches <- as.matrix(c(NA))
  # define number of control patches from percent.ctrl
  n <- round(total.patches*percent.ctrl, 0) 
  
    for (i in seq(from=min.ring, to=max.ring, by=increment)) {
      
      area.buffer <- st_buffer(priority.area, i)
      area.ring <- st_difference(area.buffer, priority.area)
      waters.subset <- st_intersection(waters, area.ring)
      control.patches.num <- nrow(as.matrix(waters.subset))
      
      
      if (is.null(control.patches.num)){
        print(paste(i," isn't far enough"))
        next
        }
      else if (control.patches.num < n){ 
        print(paste(i," isn't far enough"))
        next
        }
      else if (control.patches.num >= n) {
        control.patches <- waters.subset
        break
        } 
      
    }

  return(control.patches)
}

  #example run to show how function works
if(1==2){
 
  plot(r, col=NA, legend=NA, asp=1,
       xlab="Easting (m)",ylab="Northing (m)"
       #main="Max edge distance = 5 km"
  )
  #points(waters, pch=16)
  points(pts.data$x_coords,pts.data$y_coords,
         pch=1,col="blue1")
  trial1 <- buffered.ctrl.patches()
  trial2 <- buffered.ctrl.patches(percent.ctrl = 0.2)
  trial3 <- buffered.ctrl.patches(percent.ctrl = 0.5)
  
  center.bf <- st_buffer(st_as_sf(center.pt), radius)
  plot(center.bf, add=TRUE, lty=2)
  
  plot(trial3, pch=16, col="red", add=TRUE)
  plot(trial2, pch=16, col="black", add=TRUE)
  plot(trial1, pch=16, col="purple", add=TRUE)
  
  
  legend(x=1000, y=29000,
         c("Patches","50% closest","20% closest","5% closest"),
         cex=0.8,
         col=c("blue1","red","black","purple"),
         pch=c(1,16,16,16),
         bty="o")
  
  
}

  # for loop it

dim(simpatch[,,])
num.patches <- 170

  # choose the level of control effort

#perc.ctrl.patches <- 0.05 
#perc.ctrl.patches <- 0.2   
perc.ctrl.patches <- 0.5 

n <- round(num.patches*perc.ctrl.patches, 0)
n

for (i in 1:dim(simpatch)[3]){ # make sure level of control effort lines up!
  buffPatches   <-  buffered.ctrl.patches(waters=st_multipoint(simpatch[,,i]),
                        #percent.ctrl = 0.05)
                        #percent.ctrl = 0.2)          
                        percent.ctrl = 0.5)
  
  temp.buff.patches <- 
    rasterize(x=buffer(vect(buffPatches), 30),
              y=r,
              values=1,
              background=background)
  writeRaster(temp.buff.patches, # make sure to select appropriate path for control effort
              #filename = paste0("SimulationLandscape/ControlAreas/AreaProximity/buffTop5perc",i,".asc" ),
              #filename = paste0("SimulationLandscape/ControlAreas/AreaProximity/buffTop20perc",i,".asc" ),
              filename = paste0("SimulationLandscape/ControlAreas/AreaProximity/buffTop50perc",i,".asc" ),
              NAflag = -999,
              overwrite = TRUE 
)

  
  }


# Notes on HexSIM import ------------------------------------------------
  # before importing into hexSIM, check to make sure novalue is defined in ASC
  # potentially need to reset to -999 if 'NAflag' argument not set correctly in writeraster

  # HexSIM is creating issues with single hexcell sized options for the created Hexmaps
  # from the rasters being imported. One option that is possible is to reduce cell size on import objects
  # this really only maters on the waters files

  # one option as reducing the grain of the landscape quickly generates ASC files that are 
  # too large to open, is to set the landscape cell size to about 5 m (6000 rows, instead of the 1000 rows = 30 meters)
  # then go into the landscape and manually augment the waters post-hexsim import.
  # this will need to be done on all waters the same way, so we choose a rule that was, eliminate the top cell in dual cell waters
  # this happens because the water falls across 2 hexcells in the map in hexsim. There is no way around
  # this for efficient import based on the set rules hexsim works by.
  # We do not expect folks to try and recreate the data via this script, but they can if must. 
  # We have included all location data and files in the HexSIM directory in associated Appendices.

  # in the hexsim 

# Defining Location Trait Data --------------------------------------------

### Simulation Waters object based traits

Location_Traits<- data.frame(Name=c("Matrix",1:num.patches),
                                   Threshold=c("-Infinity",1:num.patches))
head(Location_Traits)
write.csv(Location_Traits, 
          file="SimulationLandscape/LocationTraitfilesinfo.csv",
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE)


# after this you must go into the file directly and delete the first line to get
# rid of the first row, as it won't delete it for some reason
# make sure to use notepad bc it doesn't work well in excel or spreadsheet programs

```

## Part 2: Combining Census Files

```r

# Drake et al. 2024: sIBMS for LICA control 
# Appendix S3 Pt2 
# Combining census files from HexSIM output
# For SMARTSIM project in Dr. Meryl Mims Lab
# March 2024
# Copyright J.Drake 2024

### This document is for data wrangling of the hexsim output ###

# basic strategy
# load data of relevant census
# make sure dataframe is right shape and contents
# assign individual id based on scenario, replicate, census
# combine data 
# write to csv


library("terra")
library("dplyr")
library("tidyr")
library("sf")
library("ggplot2")
library("tidyterra")
library("reshape")


## Defining which Scenarios/sets/runs
# Scenario refers to type of treatment
# Scenario Descriptions
#A - No treatment
#B - Btwn centrality 
#C - Clusters
#D - Proximity
#E - Eigenvalue centrality

# effort refers to the amount of treatment applied
# iter refers to different iteration of habitat distribution
# replicate refers to each replicate within each hab dist iteration 
  # also known as 'Run' in the first column of the data files


iter <- 1     
itermax <- 10 # different landscapes
replicate <- 1 # Run in data files
replicatemax <- 100

#censusnum <- 5  #0: census before repro (group X location) [now in for loop]
#1: census after repro (group x location) 
#2: census before removal in Priority zone (group x priority zone) 
#3: census before removal in Removal areas (group x removal area)
#4: census after removal in priority zone (group x priority zone)
#5: census after over-wintering survival (group X stage)

# we did not use all census data, but focused on a subset to answer the questions pertinent to manuscript

#### setup workspace
  
# directories # select for the dispersal scenario you are working with
  #Disp300Hex
#dir <- "Disp234Hex\\" # insert your path to output in hexsim folder
  #Disp234Hex
dir <- "Disp234Hex\\"
  #Disp167Hex
#dir <- "Disp167Hex\\"

# for each combination Run section 1 once and first to start the first replicate
# the run Section 2 to finish the replicates and write out the results

#### Section 1 ----

scenario <- "A"  # "A" "B" "C" "D" "E"  [ "F" "G" "H" "I"=sensitivity runs (234 Hex only) ]
effort <- "00"  # "00" "05" "20" "50"   [ scenarios A, F, G, H, I only use "00" ]

for(censusnum in 0:5){
# initial data
workspace <- paste0(dir,"Scenario",scenario,"\\Results\\",  # Scenario folder
                    "LICA_Scenario",scenario,"_Perc",effort,
                    "_Iter",iter,"\\",                      # Iteration folder
                    "LICA_Scenario",scenario,"_Perc",effort,
                    "_Iter",iter,"-[",replicate,"]\\",      # Replicate folder
                    "LICA_Scenario",scenario,"_Perc",effort,
                    "_Iter",iter,".",censusnum,".csv")      # Individual file

data <- read.csv(workspace)     #read in census file(s)

data <- data.frame(lapply(data, "length<-", 101))
data$Run <- 1 # this is the same as replicate, i in the loops
data$Time.Step <-0:100
data[is.na(data)] <- 0

iteration <- c(1)              #use 1 to build initial file
data <- cbind(data, iteration) #adding unique identifier for each replicate
Summary <- data                 #finish creating building intial file as "Summary" database



## First for loop is to finish replicates within set 1; 
## otherwise first replicate is counted twice if 1:replicatemax
for (i in 2:replicatemax)
  {
  workspace <- paste0(dir,"Scenario",scenario,"\\Results\\",  # Scenario folder
                      "LICA_Scenario",scenario,"_Perc",effort,
                      "_Iter",iter,"\\",                      # Iteration folder
                      "LICA_Scenario",scenario,"_Perc",effort,
                      "_Iter",iter,"-[",i,"]\\",              # Replicate folder
                      "LICA_Scenario",scenario,"_Perc",effort,
                      "_Iter",iter,".",censusnum,".csv")      # Individual file
  
  
  
  data <- read.csv(workspace)
  
  data <- data.frame(lapply(data, "length<-", 101))
  data$Run <- i # this is the same as replicate, i in the loops
  data$Time.Step <-0:100
  data[is.na(data)] <- 0
  
  iteration <- c(1)            #keep identifier at 1
  data <- cbind(data, iteration)
  Summary <- rbind(Summary, data)
}

#### Section 2 ----

## Nested for loops to cycle through setS 2-max (10 in this case) and all replicates
for (j in 2:itermax)
{
  for (i in 1:replicatemax)
  {
    workspace <- paste0(dir,"Scenario",scenario,"\\Results\\",  # Scenario folder
                        "LICA_Scenario",scenario,"_Perc",effort,
                        "_Iter",j,"\\",                      # Iteration folder
                        "LICA_Scenario",scenario,"_Perc",effort,
                        "_Iter",j,"-[",i,"]\\",              # Replicate folder
                        "LICA_Scenario",scenario,"_Perc",effort,
                        "_Iter",j,".",censusnum,".csv")      # Individual file
    
    data <- read.csv(workspace)
    
    data <- data.frame(lapply(data, "length<-", 101))
    data$Run <- i # this is the same as replicate, i in the loops
    data$Time.Step <-0:100
    data[is.na(data)] <- 0
    
    iteration <- c(j)          #set identifier to j to match iter
    data <- cbind(data, iteration)
    Summary <- rbind(Summary, data)
  }
}

## Add scenario/treatment level to the data set
Treatment <- scenario
Summary <- cbind(Summary, Treatment)

## Export for downstream analyses


results <- paste0(dir,"CompiledSimResults\\","LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

write.csv(Summary, results)
}

```

## Part 3: Down Stream Analyses

```r
# Drake et al. 2024: sIBMS for LICA control 
# Appendix S3 Pt3 
# Compiling, wrangling, and analyzing combined census files
#
# For SMARTSIM project in Dr. Meryl Mims Lab
# March 2024
# Copyright J.Drake 2024

### This document is for data wrangling of the combined hexsim output ###
### and for the downstream data analysis after data has been combined ###

library("dplyr")
library("tidyr")
library("ggplot2")
library("ggridges")

# a brief description of Census order in our HexSIM model

#Census 0: before repro (group X location)
#Census 1: census after repro (group x location) 
#Census 2: census before removal in Priority zone (group x priority zone) 
  #0:Floater, matrix
  #1:Floater, 1
  #2:Group member, matrix
  #3:Group member, 1<- this is umber in patches in priority area after removal
#Census 3: census before removal in Removal areas (group x removal area)
  #0:Floater, matrix
  #1:Floater, 1
  #2:Groupd member, matrix
  #3:Group member, 1 <- this is number in patches in the Removal Area 
#Census 4: census after removal in priority zone (group x priority zone)
  #0:Floater, matrix
  #1:Floater, 1
  #2:Groupd member, matrix
  #3:Group member, 1<- this is number in patches in priority area after removal
#Census 5: census after over-wintering survival (group X stage)



# Global Objects ----------------------------------------------------------

  ## directory
main.dir <- ""
results.dir <- #"SmartSIMWorkDirectory\\Disp300Hex\\CompiledSimResults\\"
  "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\"
  #"SmartSIMWorkDirectory\\Disp167Hex\\CompiledSimResults\\"

startgen <- 30 # keep at -1 to actual start generation i.e. if 30, starts gen =31



# Pond Occupancy ----------------------------------------------------------
# Uses Census # 0 - group X location

  ## we will first compile A X 00 and then run the for loop for the rest of the scenarios

scenario <- "A" # A for baseline  
                # F, G, H, I for sensitivity on this first for loop only       
effort <- "00"  
censusnum <-  0

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

Summaryponds <- read.csv(results)

## Convert pond counts to presence-absence data
Summarybinary <- Summaryponds[,179:349] # make database with only pond columns, and only locations for group members
                                        # CAUTION: this indexing only works for our event sequence and combination of traits
                                        # you would need to change to reflect your output
Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep

PondsOccupied

## Create new Summaryponds dataframe with Run, Time Step, and pond occupancy summed across ponds
Summaryoccupied <- data.frame(as.factor(Summaryponds$Run),                       
                              as.numeric(Summaryponds$Time.Step),
                              as.integer(PondsOccupied),
                              scenario,
                              effort)

colnames(Summaryoccupied) <- c("Run","Time.Step","PondsOccupied", "Treatment", "Effort") 


## Write to .csv - ponds occupied per timestep for all replicates for only A X 00
#occupancy results directory
if(1==2){
  write.csv(Summaryoccupied, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario",scenario,
                                    "_Perc",effort,"_Census",censusnum,
                                    "_SummaryOccupied.csv"), row.names=FALSE)
}

### after running the A X 00 code snippet above, run the  snippets/for loop below
  ## you should only have to run it once and it should spit out the right info 
  ## as long as the directories are accurate
SummaryoccupiedA <- Summaryoccupied 
scenarios <- LETTERS[2:5]
efforts <- c("05", "20", "50")

for (effort in efforts) {
  SummaryoccupiedAll <- SummaryoccupiedA 
for (i in scenarios){
  
  scenario <- i
  censusnum <- 0
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort,"_Census",censusnum,
                    "_ResultsCompiled.csv")
  
  Summaryponds <- read.csv(results)
  Summarybinary <- Summaryponds[,179:349] #make database with only pond columns, and only locations for group members
  Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
  PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep
  
  Summaryoccupied <- data.frame(as.factor(Summaryponds$Run),                       
                                as.numeric(Summaryponds$Time.Step),
                                as.integer(PondsOccupied),
                                scenario,
                                effort)
  
  colnames(Summaryoccupied) <- c("Run","Time.Step","PondsOccupied", "Treatment", "Effort") #change name to replicate maybe?
  
  SummaryoccupiedAll<-rbind(SummaryoccupiedAll, Summaryoccupied)
  
}

unique(SummaryoccupiedAll$Treatment)
SummaryoccupiedAll$Treatment <- as.factor(SummaryoccupiedAll$Treatment)
levels(SummaryoccupiedAll$Treatment) <- c("Control", "Betweenness", "Clusters", "PZ Buffer", "Eigenvalue")


if(1==2){ # turn on to write out compiled results
  write.csv(SummaryoccupiedAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                       "_Perc",effort,"_Census",censusnum,
                                       "_SummaryOccupiedAll.csv"), row.names=FALSE)
}



}


# Total Population --------------------------------------------------------
#Uses census 5  and can use Population.Size column

scenario <- "A" # A for baseline  
                # F, G, H, I for sensitivity in augmented for loop below        
effort <- "00"  
censusnum <-  5

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryPopsA <- read.csv(results)
SummaryPopsA <- data.frame(as.factor(SummaryPopsA$Run),                       
                          as.numeric(SummaryPopsA$Time.Step),
                          as.numeric(SummaryPopsA$Population.Size),
                          as.integer(SummaryPopsA$iteration),
                          scenario,
                          effort)
colnames(SummaryPopsA) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 


SummaryPopAll <- SummaryPopsA 
scenarios <- LETTERS[2:5]
efforts <- c("05", "20", "50")

for (effort in efforts) {
  SummaryPopAll <- SummaryPopsA 
for (i in scenarios){
  
  scenario <- i
  censusnum <- 0
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort,"_Census",censusnum,
                    "_ResultsCompiled.csv")
  
  SummaryPops <- read.csv(results)
  SummaryPops <- data.frame(as.factor(SummaryPops$Run),                       
                            as.numeric(SummaryPops$Time.Step),
                            as.numeric(SummaryPops$Population.Size),
                            as.integer(SummaryPops$iteration),
                            scenario,
                            effort)
  
  colnames(SummaryPops) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 
  
  SummaryPopAll<-rbind(SummaryPopAll, SummaryPops)
  
}

unique(SummaryPopAll$Treatment)
SummaryPopAll$Treatment <- as.factor(SummaryPopAll$Treatment)
levels(SummaryPopAll$Treatment) <- c("Control", "Betweenness", "Clusters", "PZ Buffer", "Eigenvalue")

if(1==2){ # turn on to write out combined data
  write.csv(SummaryPopAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                       "_Perc",effort,"_Census",censusnum,
                                       "_SummaryPopAll.csv"), row.names=FALSE)
}


}

# Priority Zone Population ------------------------------------------------
# Uses Census 4 - group x priority zone after removal

scenario <- "A"       
effort <- "00"  
censusnum <-  4

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryPZ <- read.csv(results)
SummaryPZ <- data.frame(as.factor(SummaryPZ$Run),                       
                        as.numeric(SummaryPZ$Time.Step),
                        as.numeric(SummaryPZ$Trait.Index..3),
                        as.integer(SummaryPZ$iteration),
                        scenario,
                        effort)
colnames(SummaryPZ) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 


SummaryPZA <- SummaryPZ 

scenarios <- LETTERS[2:5]
efforts <- c("05", "20", "50")

for (effort in efforts) {
  SummaryPZAll <- SummaryPZA
for (i in scenarios){
  
  scenario <- i
  censusnum <- 4
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort,"_Census",censusnum,
                    "_ResultsCompiled.csv")
  
  SummaryPZ <- read.csv(results)
  SummaryPZ <- data.frame(as.factor(SummaryPZ$Run),                       
                          as.numeric(SummaryPZ$Time.Step),
                          as.numeric(SummaryPZ$Trait.Index..3),
                          as.integer(SummaryPZ$iteration),
                          scenario,
                          effort)
  
  colnames(SummaryPZ) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 
  
  SummaryPZAll<-rbind(SummaryPZAll, SummaryPZ)
  
}

unique(SummaryPZAll$Treatment)
SummaryPZAll$Treatment <- as.factor(SummaryPZAll$Treatment)
levels(SummaryPZAll$Treatment) <- c("Control", "Betweenness", "Clusters", "PZ Buffer", "Eigenvalue")

ggplot(data=SummaryPZAll[which(SummaryPZAll$Time.Step>=31),], aes(y=Population,x=Time.Step)) +
  geom_line(aes(color=as.factor(Run)), show.legend = FALSE, alpha=0.25) +
  stat_summary(fun=mean, geom="line")+
  geom_smooth()+
  geom_vline(aes(xintercept=50), colour="black", linetype="dashed") +
  facet_wrap(as.factor(SummaryPZAll[which(SummaryPZAll$Time.Step>=31),]$Treatment), nrow=2)+
  theme_classic()

if(1==2){ # turn on to write out combined data
  write.csv(SummaryPZAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                  "_Perc",effort,"_Census",censusnum,
                                  "_SummaryPZAll.csv"), row.names=FALSE)
}
}




# Total Animals removed ---------------------------------------------------
# Uses Census 3 

scenario <- "A"       
effort <- "00"  
censusnum <-  3
results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryRA <- read.csv(results)
SummaryRA <- data.frame(as.factor(SummaryRA$Run),                       
                        as.numeric(SummaryRA$Time.Step),
                        as.numeric(SummaryRA$Trait.Index..3),
                        as.integer(SummaryRA$iteration),
                        scenario,
                        effort)
colnames(SummaryRA) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 


SummaryRAA <- SummaryRA 

scenarios <- LETTERS[2:5]
efforts <- c("05", "20", "50")

for (effort in efforts) {
  SummaryRAAll <- SummaryRAA
for (i in scenarios){
  
  scenario <- i
  censusnum <- 3
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort,"_Census",censusnum,
                    "_ResultsCompiled.csv")
  
  SummaryRA <- read.csv(results)
  SummaryRA <- data.frame(as.factor(SummaryRA$Run),                       
                          as.numeric(SummaryRA$Time.Step),
                          as.numeric(SummaryRA$Trait.Index..3),
                          as.integer(SummaryRA$iteration),
                          scenario,
                          effort)
  colnames(SummaryRA) <- c("Run","Time.Step","Population","Iteration", "Treatment", "Effort") 
  
  SummaryRAAll<-rbind(SummaryRAAll, SummaryRA)
  
}

unique(SummaryRAAll$Treatment)
SummaryRAAll$Treatment <- as.factor(SummaryRAAll$Treatment)
levels(SummaryRAAll$Treatment) <- c("Control", "Betweenness", "Clusters", "PZ Buffer", "Eigenvalue")

ggplot(data=SummaryRAAll[which(SummaryRAAll$Time.Step>=31),], aes(y=Population,x=Time.Step)) +
  geom_line(aes(color=as.factor(Run)), show.legend = FALSE, alpha=0.25) +
  stat_summary(fun=mean, geom="line")+
  geom_smooth()+
  geom_vline(aes(xintercept=50), colour="black", linetype="dashed") +
  facet_wrap(as.factor(SummaryRAAll[which(SummaryRAAll$Time.Step>=31),]$Treatment), nrow=2)+
  theme_classic()

if(1==1){
  write.csv(SummaryRAAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                 "_Perc",effort,"_Census",censusnum,
                                 "_SummaryRAAll.csv"), row.names=FALSE)
}

}





# Sensitivity Analyses ----------------------------------------------------
  # compiling data files for sensitivity analyses
  #* Total Population -----
#Uses census 5 and can use Population.Size column
if(1==2){
main.dir <- "directory"  
results.dir <-  "Disp234Hex\\CompiledSimResults\\"
}

scenario <- "A" # A for baseline  
# F, G, H, I for sensitivity in augmented for loop below        
effort <- "00"  
censusnum <-  5

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryPopsA <- read.csv(results)
SummaryPopsA <- data.frame(as.factor(SummaryPopsA$Run),                       
                           as.numeric(SummaryPopsA$Time.Step),
                           as.numeric(SummaryPopsA$Population.Size),
                           as.integer(SummaryPopsA$iteration),
                           scenario,
                           effort)
colnames(SummaryPopsA) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 


SummaryPopAll <- SummaryPopsA 

scenarios <- LETTERS[6:9] # sensitivity analyses F, G, H, & I
efforts <- "00" # for sensitivities

for (effort in efforts) {
  SummaryPopAll <- SummaryPopsA 
  for (i in scenarios){
    
    scenario <- i
    censusnum <- 5
    
    results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                      "_Perc",effort,"_Census",censusnum,
                      "_ResultsCompiled.csv")
    
    SummaryPops <- read.csv(results)
    SummaryPops <- data.frame(as.factor(SummaryPops$Run),                       
                              as.numeric(SummaryPops$Time.Step),
                              as.numeric(SummaryPops$Population.Size),
                              as.integer(SummaryPops$iteration),
                              scenario,
                              effort)
    
    colnames(SummaryPops) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 
    
    SummaryPopAll<-rbind(SummaryPopAll, SummaryPops)
    
  }
  


if(1==2){ #Turn on to write out data
  write.csv(SummaryPopAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                  "_Perc",effort,"_Census",censusnum,
                                  "_SummaryPopSensitivity.csv"), row.names=FALSE)
}

}


  #* Total Occupancy ----

if(1==2){
  main.dir <- "dir\\"  
  results.dir <-  "Disp234Hex\\CompiledSimResults\\"
}

scenario <- "A" # A for baseline  
# F, G, H, I for sensitivity in augmented for loop below        
effort <- "00"  
censusnum <-  0

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryOccsA <- read.csv(results)
SummaryOccsA <- data.frame(as.factor(SummaryPopsA$Run),                       
                           as.numeric(SummaryPopsA$Time.Step),
                           as.numeric(SummaryPopsA$Population.Size),
                           as.integer(SummaryPopsA$iteration),
                           scenario,
                           effort)
colnames(SummaryPopsA) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 


SummaryPopAll <- SummaryPopsA 

scenarios <- LETTERS[6:9] # 
efforts <- "00" # for sensitivities

for (effort in efforts) {
  SummaryPopAll <- SummaryPopsA 
  for (i in scenarios){
    
    scenario <- i
    censusnum <- 5
    
    results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                      "_Perc",effort,"_Census",censusnum,
                      "_ResultsCompiled.csv")
    
    SummaryPops <- read.csv(results)
    SummaryPops <- data.frame(as.factor(SummaryPops$Run),                       
                              as.numeric(SummaryPops$Time.Step),
                              as.numeric(SummaryPops$Population.Size),
                              as.integer(SummaryPops$iteration),
                              scenario,
                              effort)
    
    colnames(SummaryPops) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 
    
    SummaryPopAll<-rbind(SummaryPopAll, SummaryPops)
    
  }
  
  
  
  if(1==2){ #turn on to write out data
    write.csv(SummaryPopAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                    "_Perc",effort,"_Census",censusnum,
                                    "_SummaryPopSensitivity.csv"), row.names=FALSE)
  }
  
}

  #* Priority Zone Population -----
    # Uses Census 4 - group x priority zone after removal

scenario <- "A"       
effort <- "00"  
censusnum <-  4

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

SummaryPZ <- read.csv(results)
SummaryPZ <- data.frame(as.factor(SummaryPZ$Run),                       
                        as.numeric(SummaryPZ$Time.Step),
                        as.numeric(SummaryPZ$Trait.Index..3),
                        as.integer(SummaryPZ$iteration),
                        scenario,
                        effort)
colnames(SummaryPZ) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 


SummaryPZA <- SummaryPZ 

scenarios <- LETTERS[6:9] # 9 when I finishes up; for sensitivities
efforts <- "00" # for sensitivities

for (effort in efforts) {
  SummaryPZAll <- SummaryPZA
  for (i in scenarios){
    
    scenario <- i
    censusnum <- 4
    
    results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                      "_Perc",effort,"_Census",censusnum,
                      "_ResultsCompiled.csv")
    
    SummaryPZ <- read.csv(results)
    SummaryPZ <- data.frame(as.factor(SummaryPZ$Run),                       
                            as.numeric(SummaryPZ$Time.Step),
                            as.numeric(SummaryPZ$Trait.Index..3),
                            as.integer(SummaryPZ$iteration),
                            scenario,
                            effort)
    
    colnames(SummaryPZ) <- c("Run","Time.Step","Population","Iteration", "Scenario", "Effort") 
    
    SummaryPZAll<-rbind(SummaryPZAll, SummaryPZ)
    
  }
  

  if(1==2){ #Turn on to write out data
    write.csv(SummaryPZAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                   "_Perc",effort,"_Census",censusnum,
                                   "_SummaryPZSensitivities.csv"), row.names=FALSE)
  }
}


#* PZ Occupancy ----

  # this will be covered in a separate file as the processing takes a few extra steps

# Data Analysis and Statistical Tests --------------------

# Mean Population Last 20 Time Steps (Total landscape and PA) -----------------------
#Uses Census 4 for the PA pop
#Uses census 5  and can use Population.Size column

#The different data to consider
  # All pop data
  # PZ (Protected Area (formerly Priority Zone naming convention)) data
  # RA all = removal area data

#combine dat for presenting more info  
df_list <- list(SummaryPopAll=SummaryPopAll,
                SummaryPZAll=SummaryPZAll,
                SummaryRAAll=SummaryRAAll)

SummaryAll<- bind_rows(df_list, .id = "Data_set")
head(SummaryAll)

SummaryAll$Data_set <- as.factor(SummaryAll$Data_set)
SummaryAll <- SummaryAll[which(SummaryAll$Treatment!="Control"),]

# Combine data based on treatment/effort levels ---------------------------

  ###
  ## Run all dispersal levels above before running the below code !!!###
  ###

  ## The below code is to combine each of the long form data to allow statistical analysis


# census  0 and 5 -> total pop are the same
# census 0 -> ponds occupied
# census 4 -> PZ pop

 #
main.dir <- "dir\\"  
results.dir <- "SmartSIMWorkDirectory\\Disp300Hex\\CompiledSimResults\\SummarizedResults\\"
  #"SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\SummarizedResults\\"
  #"SmartSIMWorkDirectory\\Disp167Hex\\CompiledSimResults\\SummarizedResults\\"

censusnum <- c(0,2,3,4,5) #
scenario <- "All"
effort <- c("00", "05", "20", "50")
item <- c("OccupiedAll","PopAll","PZAll","RAAll")

#* Total population numbers ####

{  
results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[2],"_Census",censusnum[1],
                  "_Summary",item[2],".csv")

SummaryPop_5perc <- read.csv(results)


results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[3],"_Census",censusnum[1],
                  "_Summary",item[2],".csv")

SummaryPop_20perc <- read.csv(results)


results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[4],"_Census",censusnum[1],
                  "_Summary",item[2],".csv")

SummaryPop_50perc <- read.csv(results)

SummaryAllPop <- rbind(
  SummaryPop_5perc,
  SummaryPop_20perc[SummaryPop_20perc$Effort > 1 ,],
  SummaryPop_50perc[SummaryPop_50perc$Effort > 1 ,]
)

str(SummaryAllPop)
} 

if(1==2){
  write.csv(SummaryAllPop, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                 "_Perc","All","_Census",censusnum[1],
                                 "_PopAll.csv"), row.names=FALSE)
}

temp.data<-SummaryAllPop %>% 
    filter(Time.Step > 80) %>% 
    group_by(Treatment, Effort) %>% 
    summarise(Means = mean(Population),
              StDev = sd(Population))

temp.data[7,3] # used for % difference calcs

SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[7,3]))/((mean(Population)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[7,4]))/((sd(Population)+(temp.data[7,4]))/2)
            )

SummaryAllPopPercDiff<- SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[7,3]))/((mean(Population)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[7,4]))/((sd(Population)+(temp.data[7,4]))/2)
  )

#issues exporting tibble
SummaryAllPopPercDiff <- as.data.frame(SummaryAllPopPercDiff)
SummaryAllPopPercDiff$MeanPercDiff <- as.vector(SummaryAllPopPercDiff[1:13,5]$Means)
SummaryAllPopPercDiff$SDPercDiff <- as.vector(SummaryAllPopPercDiff[1:13,6]$StDev)

str(SummaryAllPopPercDiff)

if(1==2){ # turn on to write out data
  write.csv(SummaryAllPopPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                  "_Perc","All","_Census",censusnum[1],
                                  "_SummaryAllPopPercDiff.csv"), row.names=FALSE)
}

#* PZ population numbers ####
  # census 4

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[2],"_Census",censusnum[4],
                  "_Summary",item[3],".csv")

SummaryPZ_5perc <- read.csv(results)


results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[3],"_Census",censusnum[4],
                  "_Summary",item[3],".csv")

SummaryPZ_20perc <- read.csv(results)


results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort[4],"_Census",censusnum[4],
                  "_Summary",item[3],".csv")

SummaryPZ_50perc <- read.csv(results)

SummaryAllPZ <- rbind(
  SummaryPZ_5perc,
  SummaryPZ_20perc[SummaryPZ_20perc$Effort > 1 ,],
  SummaryPZ_50perc[SummaryPZ_50perc$Effort > 1 ,]
)

str(SummaryAllPZ)

if(1==2){
  write.csv(SummaryAllPZ, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                  "_Perc","All","_Census",censusnum[4],
                                  "_PZAll.csv"), row.names=FALSE)
}


temp.data<-SummaryAllPZ %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population))

temp.data[7,3] # used for % difference calcs

SummaryAllPZ %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[7,3]))/((mean(Population)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[7,4]))/((sd(Population)+(temp.data[7,4]))/2)
  )

SummaryAllPZPercDiff<- SummaryAllPZ %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[7,3]))/((mean(Population)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[7,4]))/((sd(Population)+(temp.data[7,4]))/2)
  )
#issues exporting tibble
SummaryAllPZPercDiff <- as.data.frame(SummaryAllPZPercDiff)
SummaryAllPZPercDiff$MeanPercDiff <- as.vector(SummaryAllPZPercDiff[1:13,5]$Means)
SummaryAllPZPercDiff$SDPercDiff <- as.vector(SummaryAllPZPercDiff[1:13,6]$StDev)

str(SummaryAllPZPercDiff)

if(1==2){
  write.csv(SummaryAllPZPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                          "_Perc","All","_Census",censusnum[1],
                                          "_SummaryAllPZPercDiff.csv"), row.names=FALSE)
}


#* Occupancy numbers ####
  # census 0 
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort[2],"_Census",censusnum[1],
                    "_Summary",item[1],".csv")
  
  SummaryOcc_5perc <- read.csv(results)
  
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort[3],"_Census",censusnum[1],
                    "_Summary",item[1],".csv")
  
  SummaryOcc_20perc <- read.csv(results)
  
  
  results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                    "_Perc",effort[4],"_Census",censusnum[1],
                    "_Summary",item[1],".csv")
  
  SummaryOcc_50perc <- read.csv(results)
  
  SummaryAllOcc <- rbind(
    SummaryOcc_5perc,
    SummaryOcc_20perc[SummaryOcc_20perc$Effort > 1 ,],
    SummaryOcc_50perc[SummaryOcc_50perc$Effort > 1 ,]
  )
  
  str(SummaryAllOcc)
 

if(1==2){
  write.csv(SummaryAllOcc, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                  "_Perc","All","_Census",censusnum[1],
                                  "_OccAll.csv"), row.names=FALSE)
}

temp.data<-SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied))

temp.data[7,3] # used for % difference calcs

SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied),
            MeanPercDiff = 100*(mean(PondsOccupied)-(temp.data[7,3]))/((mean(PondsOccupied)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(PondsOccupied)-(temp.data[7,4]))/((sd(PondsOccupied)+(temp.data[7,4]))/2)
  )

SummaryAllOccPercDiff <- SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied),
            MeanPercDiff = 100*(mean(PondsOccupied)-(temp.data[7,3]))/((mean(PondsOccupied)+(temp.data[7,3]))/2),
            SDPercDiff = 100*(sd(PondsOccupied)-(temp.data[7,4]))/((sd(PondsOccupied)+(temp.data[7,4]))/2)
  )


#issues exporting tibble
SummaryAllOccPercDiff <- as.data.frame(SummaryAllOccPercDiff)
SummaryAllOccPercDiff$MeanPercDiff <- as.vector(SummaryAllOccPercDiff[1:13,5]$Means)
SummaryAllOccPercDiff$SDPercDiff <- as.vector(SummaryAllOccPercDiff[1:13,6]$StDev)

str(SummaryAllOccPercDiff)

if(1==2){
  write.csv(SummaryAllOccPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                         "_Perc","All","_Census",censusnum[1],
                                         "_SummaryAllOccPercDiff.csv"), row.names=FALSE)
}

#* PZ Occupancy numbers ####
  # this is different to other sources based on how I had to manipulate the data
  # See PZOccResults.R before running this section

# "Output//Disp234Hex//Perc05_Census0_PZOcc.csv"

#censusnum <- c(0,2,3,4) #
#scenario <- "All"
#effort <- c("00", "05", "20", "50")
#item <- c("OccupiedAll","PopAll","PZAll","RAAll")


#starting with 234
if(1=2){
main.dir <- "dir\\"  
results.dir <- "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\SummarizedResults\\"

SummaryPZOcc_05perc <- read.csv("Output//Disp234Hex//Perc05_Census0_PZOcc.csv")
SummaryPZOcc_20perc <- read.csv("Output//Disp234Hex//Perc20_Census0_PZOcc.csv")
SummaryPZOcc_50perc <- read.csv("Output//Disp234Hex//Perc50_Census0_PZOcc.csv")
}
#300 hex disp
if(1=2){
main.dir <- "dir\\"  
results.dir <- "SmartSIMWorkDirectory\\Disp300Hex\\CompiledSimResults\\SummarizedResults\\"

SummaryPZOcc_05perc <- read.csv("Output//Disp300Hex//Perc05_Census0_PZOcc.csv")
SummaryPZOcc_20perc <- read.csv("Output//Disp300Hex//Perc20_Census0_PZOcc.csv")
SummaryPZOcc_50perc <- read.csv("Output//Disp300Hex//Perc50_Census0_PZOcc.csv")
}
#167 hex disp
if(1=2){
main.dir <- "dir\\"  
results.dir <- "SmartSIMWorkDirectory\\Disp167Hex\\CompiledSimResults\\SummarizedResults\\"

SummaryPZOcc_05perc <- read.csv("Output//Disp167Hex//Perc05_Census0_PZOcc.csv")
SummaryPZOcc_20perc <- read.csv("Output//Disp167Hex//Perc20_Census0_PZOcc.csv")
SummaryPZOcc_50perc <- read.csv("Output//Disp167Hex//Perc50_Census0_PZOcc.csv")
}

SummaryAllPZOcc <- rbind(
  SummaryPZOcc_05perc,
  SummaryPZOcc_20perc[SummaryPZOcc_20perc$Effort > 1 ,], # prevent deplication of Baseline data in each file
  SummaryPZOcc_50perc[SummaryPZOcc_50perc$Effort > 1 ,] # prevent deplication of Baseline data in each file
)

str(SummaryAllPZOcc)
 

temp.data<-SummaryAllPZOcc %>% 
  filter(X > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(PZOccRate),
            StDev = sd(PZOccRate))

temp.data # used for % difference calcs

SummaryAllPZOcc %>% 
  filter(X > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(PZOccRate),
            StDev = sd(PZOccRate),
            MeanPercDiff = 100*(mean(PZOccRate)-(temp.data[1,3]))/((mean(PZOccRate)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(PZOccRate)-(temp.data[1,4]))/((sd(PZOccRate)+(temp.data[1,4]))/2)
  )

SummaryAllPZOccPercDiff <- SummaryAllPZOcc %>% 
  filter(X > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(PZOccRate),
            StDev = sd(PZOccRate),
            MeanPercDiff = 100*(mean(PZOccRate)-(temp.data[1,3]))/((mean(PZOccRate)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(PZOccRate)-(temp.data[1,4]))/((sd(PZOccRate)+(temp.data[1,4]))/2)
  )

#issues exporting tibble
SummaryAllPZOccPercDiff <- as.data.frame(SummaryAllPZOccPercDiff)
SummaryAllPZOccPercDiff$MeanPercDiff <- as.vector(SummaryAllPZOccPercDiff[1:13,5]$Means)
SummaryAllPZOccPercDiff$SDPercDiff <- as.vector(SummaryAllPZOccPercDiff[1:13,6]$StDev)

str(SummaryAllPZOccPercDiff)

if(1==2){
  write.csv(SummaryAllPZOccPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                          "_Perc","All","_Census",censusnum[1],
                                          "_SummaryAllPZOccPercDiff.csv"), row.names=FALSE)
}




# Pop ANOVA -------------------------------------------------------------------

  # need to run the combine data section above to get to this step for each dispersal dist scenario
  # after setting directory and running the combination section at top of document, then run through ANOVAS

# we will take a look at an ANOVA via 
  #the mean value of the 20 time steps mean((time.steps 81 - 100) (all pop, pz pop, occ, pz Occ))


#* ANOVA: Mean of last 20 time step pops ####

PopAnovaData <- SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort, Run) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population))

PopAnovaData$Treatment <- as.factor(PopAnovaData$Treatment)
PopAnovaData$Effort <- as.factor(PopAnovaData$Effort)

PopAnova <- aov(formula=Means~Effort/Treatment,
                data=PopAnovaData)

summary(PopAnova)
TukeyHSD(PopAnova)


Mean20PopTukey <- TukeyHSD(PopAnova)

if(1==2){
  write.csv(Mean20PopTukey$Effort, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                      "_Perc","All","_Census",censusnum[1],
                                                      "_Mean20PopTukeyHSDmains.csv"), row.names=TRUE)
  
  write.csv(Mean20PopTukey$`Effort:Treatment`, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                       "_Perc","All","_Census",censusnum[1],
                                                       "_Mean20PopTukeyHSDinteraction.csv"), row.names=TRUE)
}

#plot(PopAnova)
plot(TukeyHSD(PopAnova), las=1)
plot(TukeyHSD(PopAnova, ordered = TRUE), las=1)


# PZ Pop ANOVA -------------------------------------------------------------------


#* ANOVA: Mean of last 20 time step pops ####

PZAnovaData <- SummaryAllPZ %>% 
  filter(Time.Step > 80) %>% 
 
  group_by(Treatment, Effort, Run) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population))

PZAnovaData$Treatment <- as.factor(PZAnovaData$Treatment)
PZAnovaData$Effort <- as.factor(PZAnovaData$Effort)

PZAnova <- aov(formula=Means~Effort/Treatment,
                data=PZAnovaData)

summary(PZAnova)
TukeyHSD(PZAnova)
TukeyHSD(PZAnova, ordered = TRUE)

Mean20PZTukey <- TukeyHSD(PZAnova)

if(1==2){
  write.csv(Mean20PZTukey$Effort, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                      "_Perc","All","_Census",censusnum[4],
                                                      "_Mean20PZTukeyHSDmains.csv"), row.names=TRUE)
  write.csv(Mean20PZTukey$`Effort:Treatment`, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                     "_Perc","All","_Census",censusnum[4],
                                                     "_Mean20PZTukeyHSDinteraction.csv"), row.names=TRUE)
}


# Occ ANOVA -------------------------------------------------------------------
   

#* ANOVA: Last 20 time step pops

SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort, Run) 

#* ANOVA: Mean of last 20 time step pops ####

OccAnovaData <- SummaryAllOcc %>% 
  filter(Time.Step > 80) %>%

  group_by(Treatment, Effort, Run) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied))

OccAnovaData$Treatment <- as.factor(OccAnovaData$Treatment)
OccAnovaData$Effort <- as.factor(OccAnovaData$Effort)


OccAnova <- aov(formula=Means~Effort/Treatment,
                data=OccAnovaData)

summary(OccAnova)
#model.tables(OccAnova, "means")
TukeyHSD(OccAnova)
TukeyHSD(OccAnova, ordered = TRUE)

Mean20OccTukey <- TukeyHSD(OccAnova)

if(1==2){
  write.csv(Mean20OccTukey$Effort, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                       "_Perc","All","_Census",censusnum[1],
                                                       "_Mean20OccTukeyHSDmains.csv"), row.names=TRUE)
  write.csv(Mean20OccTukey$`Effort:Treatment`, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                      "_Perc","All","_Census",censusnum[1],
                                                      "_Mean20OccTukeyHSDinteraction.csv"), row.names=TRUE)
  
  }





# PZ Occ ANOVA------------------------------------------------------------------

# Make sure PZ OCC dataset is setup before in the PZOccResults.R and above compilation steps
#
#starting with 234
if(1=2){
  main.dir <- "dir\\"  
  results.dir <- "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\SummarizedResults\\"
  
  SummaryPZOcc_05perc <- read.csv("Output//Disp234Hex//Perc05_Census0_PZOcc.csv")
  SummaryPZOcc_20perc <- read.csv("Output//Disp234Hex//Perc20_Census0_PZOcc.csv")
  SummaryPZOcc_50perc <- read.csv("Output//Disp234Hex//Perc50_Census0_PZOcc.csv")
}
#300 hex disp
if(1=2){
  main.dir <- "C:\\Users\\drakej\\Desktop\\"  
  results.dir <- "SmartSIMWorkDirectory\\Disp300Hex\\CompiledSimResults\\SummarizedResults\\"
  
  SummaryPZOcc_05perc <- read.csv("Output//Disp300Hex//Perc05_Census0_PZOcc.csv")
  SummaryPZOcc_20perc <- read.csv("Output//Disp300Hex//Perc20_Census0_PZOcc.csv")
  SummaryPZOcc_50perc <- read.csv("Output//Disp300Hex//Perc50_Census0_PZOcc.csv")
}
#167 hex disp
if(1=2){
  main.dir <- "C:\\Users\\drakej\\Desktop\\"  
  results.dir <- "SmartSIMWorkDirectory\\Disp167Hex\\CompiledSimResults\\SummarizedResults\\"
  
  SummaryPZOcc_05perc <- read.csv("Output//Disp167Hex//Perc05_Census0_PZOcc.csv")
  SummaryPZOcc_20perc <- read.csv("Output//Disp167Hex//Perc20_Census0_PZOcc.csv")
  SummaryPZOcc_50perc <- read.csv("Output//Disp167Hex//Perc50_Census0_PZOcc.csv")
}

SummaryAllPZOcc <- rbind(
  SummaryPZOcc_05perc,
  SummaryPZOcc_20perc[SummaryPZOcc_20perc$Effort > 1 ,], # prevent baseline valeus replication
  SummaryPZOcc_50perc[SummaryPZOcc_50perc$Effort > 1 ,] # prevent baseline valeus replication
)

str(SummaryAllPZOcc)



#* ANOVA: Mean of last 20 time Protected Area Occupancy ####

PZOccAnovaData <- SummaryAllPZOcc %>% 
  filter(X > 80) %>%
  group_by(Scenario, Effort, Replicate) %>% 
  summarise(Means = mean(PZOccRate),
            StDev = sd(PZOccRate))

PZOccAnovaData$Scenario <- as.factor(PZOccAnovaData$Scenario)
PZOccAnovaData$Effort <- as.factor(PZOccAnovaData$Effort)


PZOccAnova <- aov(formula=Means~Effort/Scenario,
                data=PZOccAnovaData)

summary(PZOccAnova)
TukeyHSD(PZOccAnova)

# mean20 refers to mean of last 20 timesteps

if(1==2){
  Mean20PZOccTukey <- TukeyHSD(PZOccAnova)
  write.csv(Mean20PZOccTukey$Effort, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                          "_Perc","All","_Census",censusnum[1],
                                          "_Mean20PZOccTukeyHSDmains.csv"), row.names=TRUE)
  write.csv(Mean20PZOccTukey$`Effort:Scenario`, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                                      "_Perc","All","_Census",censusnum[1],
                                                      "_Mean20PZOccTukeyHSDinteraction.csv"), row.names=TRUE)
  
}



# Sensitivity Scenario Stats ----------------------------------------------

####
# Getting the % difference for demographic parameter sensitivity scenarios

# set your directory
if(1=2){
main.dir <- "dir\\"  
results.dir <-  "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\SummarizedResults\\"




#* Total population numbers ####

censusnum <- 5  # 5 for pop
                # or 0 for occ
scenario <- "All"
effort <- "00"
item <- "PZAll"

{  
  results <- paste0(main.dir,results.dir,"LICA_ScenarioAll_Perc00_Census5_SummaryPopSensitivity.csv")
  
  SummaryAllPop <- read.csv(results)
  
  str(SummaryAllPop)
} 

temp.data<-SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population))

temp.data # used for % difference calcs

SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[1,3]))/((mean(Population)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[1,4]))/((sd(Population)+(temp.data[1,4]))/2)
  )

SummaryAllPopPercDiff<- SummaryAllPop %>% 
  filter(Time.Step > 80) %>% 
  group_by(Scenario, Effort) %>% 
  summarise(Means = mean(Population),
            StDev = sd(Population),
            MeanPercDiff = 100*(mean(Population)-(temp.data[1,3]))/((mean(Population)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(Population)-(temp.data[1,4]))/((sd(Population)+(temp.data[1,4]))/2)
  )

#issues exporting tibble
SummaryAllPopPercDiff <- as.data.frame(SummaryAllPopPercDiff)
SummaryAllPopPercDiff$MeanPercDiff <- as.vector(SummaryAllPopPercDiff[1:5,5]$Means)
SummaryAllPopPercDiff$SDPercDiff <- as.vector(SummaryAllPopPercDiff[1:5,6]$StDev)
SummaryAllPopPercDiff
str(SummaryAllPopPercDiff)

if(1==2){
  write.csv(SummaryAllPopPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                          "_Perc","All","_Census",censusnum,
                                          "_SummaryAllPopPercDiffSensitivity.csv"), row.names=FALSE)
}


}

#* total occupancy ----
# Uses Census # 0 - group X location

## we will first compile A X 00 and then run the for loop for the rest of the scenarios

scenario <- "A" # A for baseline  
# F, G, H, I for sensitivity on this first for loop only       
effort <- "00"  
censusnum <-  0

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

Summaryponds <- read.csv(results)

## Convert pond counts to presence-absence data
Summarybinary <- Summaryponds[,179:349] #make database with only pond columns, and only locations for group members
Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep

PondsOccupied

## Create new Summaryponds dataframe with Run, Time Step, and pond occupancy summed across ponds
Summaryoccupied <- data.frame(as.factor(Summaryponds$Run),                       
                              as.numeric(Summaryponds$Time.Step),
                              as.integer(PondsOccupied),
                              scenario,
                              effort)

colnames(Summaryoccupied) <- c("Run","Time.Step","PondsOccupied", "Treatment", "Effort") #change name to replicate maybe?

## Write to .csv - ponds occupied per timestep for all replicates for only A X 00
#occupancy results directory
if(1==2){
  write.csv(Summaryoccupied, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario",scenario,
                                    "_Perc",effort,"_Census",censusnum,
                                    "_SummarySensitivityOccupied.csv"), row.names=FALSE)
}

### after running the A X 00 code snippet above, run the  snippets/for loop below
## you should only have to run it once and it should spit out the right info 
## as long as the directories are accurate
SummaryoccupiedA <- Summaryoccupied 
scenarios <- LETTERS[6:9]
efforts <- "00"

for (effort in efforts) {
  SummaryoccupiedAll <- SummaryoccupiedA 
  for (i in scenarios){
    
    scenario <- i
    censusnum <- 0
    
    results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                      "_Perc",effort,"_Census",censusnum,
                      "_ResultsCompiled.csv")
    
    Summaryponds <- read.csv(results)
    Summarybinary <- Summaryponds[,179:349] #make database with only pond columns, and only locations for group members
    Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
    PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep
    
    Summaryoccupied <- data.frame(as.factor(Summaryponds$Run),                       
                                  as.numeric(Summaryponds$Time.Step),
                                  as.integer(PondsOccupied),
                                  scenario,
                                  effort)
    
    colnames(Summaryoccupied) <- c("Run","Time.Step","PondsOccupied", "Treatment", "Effort") #change name to replicate maybe?
    
    SummaryoccupiedAll<-rbind(SummaryoccupiedAll, Summaryoccupied)
    
  }
  
  unique(SummaryoccupiedAll$Treatment)
  SummaryoccupiedAll$Treatment <- as.factor(SummaryoccupiedAll$Treatment)
  levels(SummaryoccupiedAll$Treatment) <- c("Control", "Betweenness", "Clusters", "PZ Buffer", "Eigenvalue")
  
  
  if(1==2){
    write.csv(SummaryoccupiedAll, paste0(main.dir,results.dir,"SummarizedResults\\","LICA_Scenario","All",
                                         "_Perc",effort,"_Census",censusnum,
                                         "_SummarySensitivityOccupiedAll.csv"), row.names=FALSE)
  }
}

SummaryAllOcc <- SummaryoccupiedAll
temp.data<-SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied))

temp.data # used for % difference calcs

SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied),
            MeanPercDiff = 100*(mean(PondsOccupied)-(temp.data[1,3]))/((mean(PondsOccupied)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(PondsOccupied)-(temp.data[1,4]))/((sd(PondsOccupied)+(temp.data[1,4]))/2)
  )

SummaryAllOccPercDiff <- SummaryAllOcc %>% 
  filter(Time.Step > 80) %>% 
  group_by(Treatment, Effort) %>% 
  summarise(Means = mean(PondsOccupied),
            StDev = sd(PondsOccupied),
            MeanPercDiff = 100*(mean(PondsOccupied)-(temp.data[1,3]))/((mean(PondsOccupied)+(temp.data[1,3]))/2),
            SDPercDiff = 100*(sd(PondsOccupied)-(temp.data[1,4]))/((sd(PondsOccupied)+(temp.data[1,4]))/2)
  )


#issues exporting tibble
SummaryAllOccPercDiff <- as.data.frame(SummaryAllOccPercDiff)
SummaryAllOccPercDiff$MeanPercDiff <- as.vector(SummaryAllOccPercDiff[1:5,5]$Means)
SummaryAllOccPercDiff$SDPercDiff <- as.vector(SummaryAllOccPercDiff[1:5,6]$StDev)

str(SummaryAllOccPercDiff)

if(1==2){
  main.dir <- "C:\\Users\\drakej\\Desktop\\"  
  results.dir <-  "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\SummarizedResults\\"
  
  write.csv(SummaryAllOccPercDiff, paste0(main.dir,results.dir,"LICA_Scenario","All",
                                          "_Perc","All","_Census",censusnum[1],
                                          "_SummaryAllOccSensitivityPercDiff.csv"), row.names=FALSE)
}

```

## Part 3.5: Protected Area additional data

```r
 Drake et al. 2024: sIBMS for LICA control 
# Appendix S3 Pt3.5 
# This could be run before Pt3 or just before the Data Analysis and Statistical Tests section of Pt3
# Compiling and wrangling data to determine each PZ occ as the indexing has to be specific to each iteration
# For SMARTSIM project in Dr. Meryl Mims Lab
# March 2024
# Copyright J.Drake 2024

### This document is for data wrangling of the combined hexsim output ###

library("dplyr")
library("tidyr")

# set directory
main.dir <- "dir\\" 
# sub directy subset by dispersal
results.dir <- #"SmartSIMWorkDirectory\\Disp300Hex\\CompiledSimResults\\"
    "SmartSIMWorkDirectory\\Disp234Hex\\CompiledSimResults\\"
  #"SmartSIMWorkDirectory\\Disp167Hex\\CompiledSimResults\\"

# Pond Occupancy ----------------------------------------------------------
# Uses Census # 0 - group X location
# reworked to be able to identify each location

scenario <- "A"       
effort <- "00"  
censusnum <-  0

results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                  "_Perc",effort,"_Census",censusnum,
                  "_ResultsCompiled.csv")

Summaryponds <- read.csv(results)

## Convert pond counts to presence-absence data
Summarybinary <- Summaryponds[,180:349] #make database with only pond columns, and only locations for group members
Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
colnames(Summarybinary) <- 0:170
Summarybinary
PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep


## Create new Summaryponds dataframe with Run, Time Step, and pond occupancy summed across ponds
Summaryoccupied <- data.frame(Scenario=scenario,
                              Effort=effort,
                              Iteration=as.factor(Summaryponds$iteration),
                              Replicate=as.factor(Summaryponds$Run),                       
                              TimeStep=as.numeric(Summaryponds$Time.Step),
                              TotalOcc = as.integer(PondsOccupied)
)

Summaryoccupied <- cbind(Summaryoccupied, Summarybinary)

str(Summaryoccupied)
dim(Summaryoccupied)
#colnames(Summaryoccupied) <- c("Run","Time.Step","PondsOccupied", "Treatment", "Effort") #change name to replicate maybe?


### after running the A X 00 code snippet above, run the  snippets/for loop below
## you should only have to run it once and it should spit out the right info 
## as long as the directories are accurate
SummaryoccupiedA <- Summaryoccupied 
scenarios <- LETTERS[2:5]
efforts <- c("05", "20", "50")

for (effort in efforts) {
  SummaryoccupiedAll <- SummaryoccupiedA 
  for (i in scenarios){
    
    scenario <- i
    censusnum <- 0
    
    results <- paste0(main.dir,results.dir,"LICA_Scenario",scenario,
                      "_Perc",effort,"_Census",censusnum,
                      "_ResultsCompiled.csv")
    
    Summaryponds <- read.csv(results)
    Summarybinary <- Summaryponds[,180:349] #make database with only pond columns, and only locations for group members [179-349 in others?]
    Summarybinary[Summarybinary > 0] <- 1 #convert to presence-absence
    colnames(Summarybinary) <- 1:170
    PondsOccupied <- rowSums(Summarybinary, na.rm=TRUE) #number of ponds occupied for each timestep
    
    
    
    ## Create new Summaryponds dataframe with Run, Time Step, and pond occupancy summed across ponds
    Summaryoccupied <- data.frame(Scenario=scenario,
                                  Effort=effort,
                                  Iteration=as.factor(Summaryponds$iteration),
                                  Replicate=as.factor(Summaryponds$Run),                       
                                  TimeStep=as.numeric(Summaryponds$Time.Step),
                                  TotalOcc = as.integer(PondsOccupied)
    )
    
    Summaryoccupied <- cbind(Summaryoccupied, Summarybinary)
    
    SummaryoccupiedAll<-rbind(SummaryoccupiedAll, Summaryoccupied)
    
  }
  
  
  if(1==2){
    # change directory from Disp234Hex to Disp300Hex to Disp167Hex as needed
    write.csv(SummaryoccupiedAll, paste0("Output\\Disp167Hex\\","LICA_Scenario","All",
                                         "_Perc",effort,"_Census",censusnum,
                                         "_SummaryOccupancyAll.csv"), row.names=FALSE)
  }
  
  
  
  
  
  
  
  
}



# Subsetting to inside PZ occupancy ---------------------------------------

library("terra")
library("dplyr")
library("tidyr")
library("sf")
library("ggplot2")
library("tidyterra")

# this is a map of the different habitat prioritization methods in an example patch network
{
  set.seed(42) # set seed for reproducibility 
  
  # this is an adaptation of the landscape gen process to id those patches inside the Protected Area
  
  #* Create uniform landscape raster ----
  # Rectangle extent = 30km^2 
  # Cell size = 30 meter
  # Cell value = 1 ( or 0? )
  
  # for 30mX30m grid nrow&ncol=1000
  # for 5mX5m grid = 6000 rows
  # fpr 1x1m grid = 30,000 rows cols # super computationally intensive
  
  n.col <- 1000
  n.row <- 1000
  x.min <- 0
  y.min <- 0
  x.max <- 30000
  y.max <- 30000
  vals <- 1
  background <- -999 # for mode to import waters correctly, we will use
  # a special background for the water patch generation
  # specifically, 999
  
  r <- rast(ncol=n.col,
            nrow=n.row,
            xmin=x.min,
            xmax=x.max,
            ymin=y.min,
            ymax=y.max,
            #crs=,
            #extent=,
            #resolution=,
            vals=vals,
  )
  
  #* Create management area in center ----

  center.pt <- matrix(data=NA, nrow=1, ncol=2) 
  center.pt[,1] <- x.max/2
  center.pt[,2] <- y.max/2
  center.pt <- vect(center.pt)
  
  plot(r, asp=1)
  plot(center.pt, add=TRUE)
  radius <- sqrt((0.15*x.max*y.max)/pi) # radius set to encapsulate 15% of the total landscape area
  
  center.bf <- buffer(center.pt, radius)
  MgmtArea <- mask(x=r,
                   mask=center.bf
  )
  
  #* Create static waters ----
  # The settings for all iterations
  
  num.patches <- 170
  min.border <- 100
  max.border <- 29900
  n.iter <- 100
  
  # bucket for generated data
  
  simpatch <- array(NA, dim=c(num.patches, 2, n.iter)) # 100 patches (rows), 2 fields xy, and 100 iterations
  
  # for loop for generating landscapes 
  
  for(i in 1:n.iter) {
    simpatch[,1,i] <- round(runif(num.patches, min=min.border, max=max.border), 0)
    simpatch[,2,i] <- round(runif(num.patches, min=min.border, max=max.border), 0)
  }
  simpatch
  
  
  #* making points spatial ----
  
  dim(simpatch[,,])
  PZ.pts<- data.frame(matrix(ncol=10,nrow=170))
  
  for (i in 1:10){
    pts <- simpatch[,,i]
    pts <- as.data.frame(pts)
    intersecting.pts <- st_intersects(st_as_sf(center.bf), st_as_sf(pts, coords = c('V1', 'V2')), sparse = FALSE)
    sum(intersecting.pts)
    PZ.pts[,i]<- as.vector(intersecting.pts)
  }
  PZ.pts <- apply(as.matrix(PZ.pts),2, as.integer)
  colSums(PZ.pts)
  
  write.csv(PZ.pts, "Output//PriorityZonePts.csv")
}

# first need to subset by landscape iteration   

# for loop ready 
#PZin <- pts[,i]
#df[df$iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]

#* Disp234Occ 05 Effort ----


Disp234Occ05<- read.csv("Output//Disp234Hex//LICA_ScenarioAll_Perc05_Census0_SummaryOccupancyAll.csv")
Disp234Occ20<- read.csv("Output//Disp234Hex//LICA_ScenarioAll_Perc20_Census0_SummaryOccupancyAll.csv")
Disp234Occ50<- read.csv("Output//Disp234Hex//LICA_ScenarioAll_Perc50_Census0_SummaryOccupancyAll.csv")

Occ.DF <- Disp234Occ05

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp234Occ05[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp234Occ05[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  #(ncol(temp.occ)-6)
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  # and if you wanted to you could augment this code to spit out a different file for 
  # each landscape iteration and keep specific occupancy info for each pondXtimestep
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp234Hex//Perc05_Census0_PZOcc.csv") 


#* Disp234Occ 20% Effort ----

Occ.DF <- Disp234Occ20

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp234Occ20[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp234Occ20[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp234Hex//Perc20_Census0_PZOcc.csv")

#* Disp234Occ 50% Effort ----

Occ.DF <- Disp234Occ50

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp234Occ50[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp234Occ50[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
   PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp234Hex//Perc50_Census0_PZOcc.csv")

#* Import disp300hex data ---

Disp300Occ05<- read.csv("Output//Disp300Hex//LICA_ScenarioAll_Perc05_Census0_SummaryOccupancyAll.csv")
Disp300Occ20<- read.csv("Output//Disp300Hex//LICA_ScenarioAll_Perc20_Census0_SummaryOccupancyAll.csv")
Disp300Occ50<- read.csv("Output//Disp300Hex//LICA_ScenarioAll_Perc50_Census0_SummaryOccupancyAll.csv")


#* Disp300Occ 05 Effort ----

Occ.DF <- Disp300Occ05

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp300Occ05[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp300Occ05[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp300Hex//Perc05_Census0_PZOcc.csv") 


#* Disp300Occ 20% Effort ----

Occ.DF <- Disp300Occ20

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp300Occ20[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp300Occ20[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp300Hex//Perc20_Census0_PZOcc.csv")


#* Disp300Occ 50% Effort ----

Occ.DF <- Disp300Occ50

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp300Occ50[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp300Occ50[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp300Hex//Perc50_Census0_PZOcc.csv")


#* Import disp167hex data ---

Disp167Occ05<- read.csv("Output//Disp167Hex//LICA_ScenarioAll_Perc05_Census0_SummaryOccupancyAll.csv")
Disp167Occ20<- read.csv("Output//Disp167Hex//LICA_ScenarioAll_Perc20_Census0_SummaryOccupancyAll.csv")
Disp167Occ50<- read.csv("Output//Disp167Hex//LICA_ScenarioAll_Perc50_Census0_SummaryOccupancyAll.csv")


#* Disp167Occ 05 Effort ----

Occ.DF <- Disp167Occ05

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp167Occ05[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp167Occ05[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin) 
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp167Hex//Perc05_Census0_PZOcc.csv") 


#* Disp167Occ 20% Effort ----

Occ.DF <- Disp167Occ20

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp167Occ20[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp167Occ20[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp167Hex//Perc20_Census0_PZOcc.csv")


#* Disp167Occ 50% Effort ----

Occ.DF <- Disp167Occ50

PZOccAll <- matrix(ncol=9,nrow=0)
colnames(PZOccAll) <- c(colnames(Disp167Occ50[,1:6]),"OccinPZ", "PondsinPZ", "PZPccRate")
for (i in 1:10) {
  PZin <- PZ.pts[,i]
  
  temp.occ <- Occ.DF[Occ.DF$Iteration == i, as.logical(c(1,1,1,1,1,1,PZin))]
  OccinPZ <- rowSums(temp.occ[,7:ncol(temp.occ)])
  
  PZOcc <- Disp167Occ50[Occ.DF$Iteration == i,1:6] 
  PZOcc$OccinPZ <- OccinPZ
  PZOcc$PondsinPZ <- sum(PZin)  
  PZOcc$PZOccRate <- OccinPZ/sum(PZin)
  PZOccAll <- rbind(PZOccAll, PZOcc)
  
}

write.csv(PZOccAll, "Output//Disp167Hex//Perc50_Census0_PZOcc.csv")


```










