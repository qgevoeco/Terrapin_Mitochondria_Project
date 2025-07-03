rm(list = ls())
setwd("~/Documents/Turtles/TurtleGenomics/terpMitoHaplo_MS2022-2023/Figures")


library(leaflet)
library(leaflet.extras)
library(mapview)  #<-- requires installation of `webshot` package
# requires PhantomJS: webshot::install_phantomjs()
################################################################################

#################
# FIGURE 1 Map
#################

## Read in Wolak research group Whole Genome Sequence Data
studSites <- read.csv("wgsSampleLocations.csv", header = TRUE)



# Map
lfmp <- leaflet(data = studSites) %>%
  # background map tile
  addProviderTiles(providers$CartoDB.Voyager) %>%
#  addProviderTiles(providers$OpenStreetMap.HOT) %>% 

  # bound the map region
  fitBounds(lng1 = -88.700, lat1 = 30.44,
    lng2 = -88.01, lat2 = 30.42) %>%

  # include points where obtained samples
  addCircleMarkers(~ Longitude, ~ Latitude,
    radius = 7, color = "black", weight = 2, opacity = 1,
    fillColor = ~ clr, fillOpacity = 1,
    label = ~ Location) %>%

  # "hand draw" scale bar (from `leaflet.extras` package)
  ## use path measurement to verify line drawn is correct measure
  ## guess and check with coordinates to place and make length of line desired
  enableMeasurePath() %>%
  addPolylines(lng = c(-88.05, -88.05), lat = c(30.25, 30.322),
    weight = 8, color = "black", opacity = 1, fillOpacity = 1) %>%
  addMeasurePathToolbar(measurePathOptions(showOnHover = FALSE,
    showArea = FALSE, imperial = FALSE)) %>%
  ## label the line with the distance (use measure path tool)  
  addLabelOnlyMarkers(lng = -88.04, lat = 30.285, label = "8km", icon = FALSE,
    labelOptions = list(permanent = TRUE, direction = "auto", offset = c(31, 0),
      textsize = "16px", textOnly = TRUE)) %>%
      
  # include inset with larger geographical area        
  addMiniMap(position = "topleft", width = 250, height = 250,
    zoomLevelOffset = -7, aimingRectOptions = list(weight = 4)) %>%

  # add a scale bar, but this doesn't end up getting saved to the file    
  addScaleBar(position = "topleft",
    options = scaleBarOptions(metric = TRUE, imperial = FALSE,
      updateWhenIdle = TRUE))

# save map to file
mapshot(lfmp, file = "Fig1_wgsSampleMap.png")

################################################################################


#################
# FIGURE 2 Map
#################

## Read in the Parham et al. 2008 ND3 to ND4 data
P08 <- read.csv("Parham2008Locations.csv", header = TRUE)
  ### Remove duplicated sites
  dupes <- which(duplicated(P08$Latitude) & duplicated(P08$Longitude))
  P08 <- P08[-dupes, ]
  
### Add colors from haplotype network
ylw <- "#f1b500"
grn <- "#00ae9a"

#### Assign haplotypes:
##### A and C to "Gulf Coast" (green)
##### B, D through O to "Atlantic Coast" (yellow)
P08$HapGrp <- "AtlanticCoast"
  P08$HapGrp[which(P08$Haplotype == "A")] <- "GulfCoast"
  P08$HapGrp[which(P08$Haplotype == "C")] <- "GulfCoast"
P08$clr <- ylw
  P08$clr[which(P08$HapGrp == "GulfCoast")] <- grn

## Add together Alabama data with Parham et al. 2008
studSites2 <- rbind(studSites[, c("Latitude", "Longitude", "clr", "Location")],
  P08[, c("Latitude", "Longitude", "clr", "Location")])

## Determine map region bounds
apply(studSites2[, c("Latitude", "Longitude")], MARGIN = 2, FUN = range)


# Map
lfmp2 <- leaflet(data = studSites2) %>%
  # background map tile
  addProviderTiles(providers$CartoDB.Voyager) %>%
#  addProviderTiles(providers$OpenStreetMap.HOT) %>% 

  # bound the map region  
  fitBounds(lng1 = -98, lat1 = 39.5,
    lng2 = -64.0, lat2 = 24.0) %>%

  # include points where obtained samples
  addCircleMarkers(~ Longitude, ~ Latitude,
    radius = 7, color = "black", weight = 2, opacity = 1,
    fillColor = ~ clr, fillOpacity = 1,
    label = ~ Location) %>%
      
  # add a scale bar, but this doesn't end up getting saved to the file    
  addScaleBar(position = "topleft",
    options = scaleBarOptions(metric = TRUE, imperial = FALSE,
      updateWhenIdle = TRUE)) %>%

  # "hand draw" scale bar (from `leaflet.extras` package)
  ## use path measurement to verify line drawn is correct measure
  ## guess and check with coordinates to place and make length of line desired
  enableMeasurePath() %>%
  addPolylines(lng = rep(-75, 2), lat = c(28.1, 29.9),
    weight = 8, color = "black", opacity = 1, fillOpacity = 1) %>%
  addMeasurePathToolbar(measurePathOptions(showOnHover = FALSE,
    showArea = FALSE, imperial = FALSE)) %>%
  ## label the line with the distance (use measure path tool)  
  addLabelOnlyMarkers(lng = -75, lat = 29.0, label = "200km", icon = FALSE,
    labelOptions = list(permanent = TRUE, direction = "auto", offset = c(60, 0),
      textsize = "16px", textOnly = TRUE))

# save map to file
mapshot(lfmp2, file = "Fig2_rangeHapMap.png")

################################################################################
# Versions with pasted output
## last updated 2025 July 3
R.Version()[c("platform", "svn rev", "version.string")]
#$platform
#[1] "x86_64-pc-linux-gnu"

#$`svn rev`
#[1] "88306"

#$version.string
#[1] "R version 4.5.1 (2025-06-13)"

packageVersion("leaflet")
#[1] ‘2.2.2’
packageVersion("leaflet.extras")
#[1] ‘2.0.1’
packageVersion("webshot")
#[1] ‘0.5.5’
packageVersion("mapview")
#[1] ‘2.11.2’




