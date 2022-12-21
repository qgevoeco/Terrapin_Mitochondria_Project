###################################################################################################
# Generating Population Genetic summaries for Terrapin mitochondrial DNA
###################################################################################################
library(tidyverse)
library(ape)
library(pegas)
library(PopGenome)
library(readr)

############################################################
#Read in data and add population designations
############################################################
setwd("/Users/samweaver/Docs/TerrapinRProject/")

TerrapinGenome <- readData("Data/VCF_Directory/", format = "VCF", gffpath="Data/Mito_Annotation", include.unknown = TRUE, FAST = TRUE)
#Read in our genetic data
get.sum.data(TerrapinGenome)
Terrapin_info <- read_delim("~/Docs/TerrapinRProject/Data/PopGenomePops.txt", delim = "\t")
#Split the data into lists containing the indivs in each population
populations <- split(Terrapin_info$ind, Terrapin_info$pop)
#Add the populations to the GENOME object
Terrapins <- set.populations(TerrapinGenome, populations, diploid = F)
#Verify that the populations look correct
Terrapins@populations

#Get a bunch of statistics from the genome
Terrapin_sw_diversity <- diversity.stats(TerrapinGenome, pi = TRUE)
Terrapin_sw_neutrality <- neutrality.stats(TerrapinGenome)
Terrapin_sw_detailed_stats <- detail.stats(TerrapinGenome)

#This gives us the average number of pairwise differences in each window
Terrapin_Pi <- Terrapin_sw_diversity@nuc.diversity.within/15615
Terrapin_Tajima <- Terrapin_sw_neutrality@Tajima.D
Terrapin_theta<-Terrapin_sw_neutrality@theta_Watterson
TerrapinGenome <- neutrality.stats(TerrapinGenome) 

############################################################################################################
#Simulate population genotypes under neutral evolutionary history and compare observed/simulated Tajima's D
############################################################################################################
setwd("/Users/samweaver/Downloads/msdir")
Simulation <- MS(TerrapinGenome, niter=5000, thetaID="Tajima", neutrality=TRUE)
Stats<-MS_getStats(Simulation, locus=1, population=1)
Stats<-as.data.frame(Stats)
Tajima_quantile<-quantile(Stats$Tajima.D, c(.025, .50, .975))
hist(Stats$Tajima.D, breaks=100)

#Let's see if our values of Tajima's D lies outside the 2.5-97.5 quantile range
Tajima_quantile
Terrapin_Tajima

############################################################################################################
#Conduct window analysis- Step 1) Set up the windows
############################################################################################################
# set chromosome size
Mito_Length <- 15615
# set window size (200 bp) and window jump (50) for overlapping windows
window_size <- 200
window_jump <- 50
# use seq to find the start points of each window
window_start <- seq(from = 1, to = Mito_Length, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size
# no windows start before the end of chromosome 8
sum(window_start > Mito_Length)
# but some window stop positions do occur past the final point
sum(window_stop > Mito_Length)
# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < Mito_Length)]
window_stop <- window_stop[which(window_stop < Mito_Length)]
# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

############################################################################################################
# Step 2) Get all the stats in sliding windows when treated as one populations
############################################################################################################
Terrapin_sw_OnePop <- sliding.window.transform(TerrapinGenome, width = 200, jump = 50, type = 2)
#Get pi
Terrapin_sw_OnePop <- diversity.stats(Terrapin_sw_OnePop, pi = TRUE)
#Get F_ST
Terrapin_sw_OnePop <- F_ST.stats(Terrapin_sw_OnePop, mode = "nucleotide")
#Get Tajima's D
Terrapin_sw_OnePop <- neutrality.stats(Terrapin_sw_OnePop)
#Extract stats for visualization
nd_OnePop <- Terrapin_sw@nuc.diversity.within/200
Taj_D_OnePop<-Terrapin_sw@Tajima.D

############################################################################################################
# Step 3) Get stats separating the three main potential populations from one another
############################################################################################################
Terrapin_sw <- sliding.window.transform(Terrapins, width = 200, jump = 50, type = 2)

#Get pi (and Dxy)
Terrapin_sw <- diversity.stats(Terrapin_sw, pi = TRUE)
#Get F_ST
Terrapin_sw <- F_ST.stats(Terrapin_sw, mode = "nucleotide")
#Extract stats for visualization
nd <- Terrapin_sw@nuc.diversity.within/200

#make population name vector
pops <- c("CedarPoint", "Dauphin", "Mississippi")
# set population names
colnames(nd) <- paste0(pops, "_pi")

#extract fst values
fst <- t(Terrapin_sw@nuc.F_ST.pairwise)
fst[is.na(fst)] <- 0
fst[fst < 0] <- 0
is.matrix(fst)
mean(fst[,1])

dxy <- get.diversity(Terrapin_sw, between = T)[[2]]/200

# get column names 
x <- colnames(fst)
# does the same thing as above but by indexing the pops vector
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)
# replace forward slash
x <- sub("/", "_", x)
# look at x to confirm the replacement has occurred
x
#Change column names
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")
#Combine everything
Terp_data <- data.frame(windows, nd, fst, dxy)

#Check out nucleotide diversity stats
CP_pi<-mean(Terp_data$CedarPoint_pi)
DI_pi<-mean(Terp_data$Dauphin_pi)
MS_pi<-mean(Terp_data$Mississippi_pi)

#Check out Pairwise Dxy estimates
dxy1<-mean(Terp_data$CedarPoint_Dauphin_dxy)
dxy2<-mean(Terp_data$CedarPoint_Mississippi_dxy)
dxy3<-mean(Terp_data$Dauphin_Mississippi_dxy)

(dxy1+dxy2+dxy3)/3
# select nucleotide diversity data and calculate means
Terp_data %>% select(contains("_pi")) %>% summarise_all(mean)
# gather the data
pi_g <- Terp_data %>% select(contains("_pi")) %>% gather(key = "species", value = "pi")

#################################################################
#Visualize these patterns along the genome
#################################################################
# select data of interest
hs <- Terp_data %>% select(mid, CedarPoint_pi, Dauphin_pi, CedarPoint_Dauphin_fst, CedarPoint_Dauphin_dxy)
hs_g <- gather(hs, -mid, key = "stat", value = "value")
hs_g <- hs_g[!duplicated(hs_g), ]
# construct a plot with facets
a <- ggplot(hs_g, aes(mid, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")


















