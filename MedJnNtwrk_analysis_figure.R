rm(list = ls())
library(ape)
library(pegas)

################################################################################
################################################################################

################################
# AL and MS + rangewide ND3-ND4
################################
## Import fasta file, create haplotypes, and construct the median joining network
nd3nd4_seq <- read.FASTA("ND3_ND4_Alignment.fasta")
  # create columns for sample identification and grouping
  nms <- names(nd3nd4_seq)
  spltNms <- strsplit(nms, split = "_")
  # Following is messy since not standard formatting of text
  smplState <- sapply(spltNms, FUN = "[", i = 1)
    smplState[c(20:43, 45:47)] <- sapply(spltNms[c(20:43, 45:47)],
      FUN = "[", i = 2)
  smplLoc <- sapply(spltNms, FUN = "[", i = 2)
    smplLoc[grepl("Dauphin Island", smplLoc)] <- "Dauphin Island"
    smplLoc[9] <- "Cedar Point"
    smplLoc[c(10, 19)] <- "MS"
    smplLoc[44] <- "MD"
  smplID <- sapply(spltNms, FUN = "[", i = 3)
    smplID[c(9:10, 19)] <- c("22", "1", "2")
    smplID[c(20:47)] <- sapply(spltNms[c(20:47)], FUN = "[", i = 1) 
    smplID[44] <- "MD.1"
  smplRegion <- rep(NA, length(nms))
    smplRegion[which(smplState == "Alabama" | smplState == "Mississippi" |
      smplState == "Texas" | smplState == "Louisiana")] <- "Gulf"
    smplRegion[which(smplState == "Florida" & smplID == "5")] <- "Gulf"
    smplRegion[c(25:29)] <- "FL Keys"
    smplRegion[c(30:47)] <- "Atlantic"
    
nd3nd4_haps <- haplotype(nd3nd4_seq)
  nd3nd4_sz <- summary(nd3nd4_haps)
nd3nd4_ntwrk <- mjn(nd3nd4_haps)
  nd3nd4_ntwrkLbls <- attr(nd3nd4_ntwrk, "labels")
# re-order frequencies from haplotypes to be in order of network
nd3nd4_sz <- nd3nd4_sz[head(nd3nd4_ntwrkLbls, -1)]  #<-- head() to get rid of NA

# Get frequencies by region
nd3nd4_hapIndx <- attr(nd3nd4_haps, "index")
  names(nd3nd4_hapIndx) <- rownames(nd3nd4_haps)
nd3nd4_hapDF <- data.frame(hap = rep(names(nd3nd4_hapIndx),
                                         lengths(nd3nd4_hapIndx)),
  index = unlist(nd3nd4_hapIndx),
  row.names = NULL)
nd3nd4_hapDF$Region <- smplRegion[nd3nd4_hapDF$index] 
nd3nd4_hapDF$Loc <- smplLoc[nd3nd4_hapDF$index] 
nd3nd4_hapDF$smplID <- smplID[nd3nd4_hapDF$index] 

  
### Add colors from haplotype network
pnkCP <- "#f09abf"
fuschDI <- "#cb0051"
purMS <- "#5f0a48"
ylwAt <- "#f1b500"
grnGu <- "#00ae9a"
cyanKys <- "#00A08A"

nd3nd4_hapDF$clr <- ylwAt
  nd3nd4_hapDF$clr[which(nd3nd4_hapDF$Region == "FL Keys")] <- grnGu #cyanKys  #<-- do this to distinguish the Keys group
  nd3nd4_hapDF$clr[which(nd3nd4_hapDF$Region == "Gulf")] <- grnGu
  nd3nd4_hapDF$clr[which(nd3nd4_hapDF$Loc == "Dauphin Island")] <- fuschDI
  nd3nd4_hapDF$clr[which(nd3nd4_hapDF$Loc == "Cedar Point")] <- pnkCP
  nd3nd4_hapDF$clr[which(nd3nd4_hapDF$Loc == "MS")] <- purMS
# Make sample groupings
nd3nd4_hapDF$smplGrp <- nd3nd4_hapDF$Loc
  nd3nd4_hapDF$smplGrp[which(nd3nd4_hapDF$clr == ylwAt)] <- "Atlantic"
  nd3nd4_hapDF$smplGrp[which(nd3nd4_hapDF$clr == grnGu)] <- "Gulf"

regFreq <- with(nd3nd4_hapDF, table(hap, smplGrp))[head(nd3nd4_ntwrkLbls, -1), ]
  # convert into matrix from table
  regFreq <- matrix(regFreq, nrow = nrow(regFreq), ncol = ncol(regFreq),
    #dimnames = dimnames(regFreq))
    dimnames = list(dimnames(regFreq)[[1L]],  #<-- do this way to drop names
                    dimnames(regFreq)[[2L]]))

xyFxd <- list(x = c(-19.2703041615388, -26.2555041836221, -19.2703041615388, 
-20.3660218120617, -29.8733119223559, -42.9907370049949, -39.1679516111867, 
2.96348234812124, -0.145748991824565, 9.08139004574056, 20.8872865431029, 
19.0164083816373, 12.7794371162552, 14.6688238632113, 26.0083734559548, 
3.23782628988116, -8.72402177525614, -10.0221308951818), y = c(-4.70161373214612, 
7.71341576437304, 4.28698711594383, 11.8823009262514, 0.452980282393151, 
4.19679095716563, -11.2125713481116, 1.7348868406692, -5.31136007924131, 
-4.48467022403503, 14.7292642773873, 2.22492682865426, -10.5151846452705, 
6.767920354891, 19.8595162252866, 7.68301735747678, 9.27821515344515, 
1.82639654092778))


pegas:::plot.haploNet(nd3nd4_ntwrk,
  size = sqrt(nd3nd4_sz),
  xy = xyFxd,
  pie = regFreq,
  lwd = 2, 
  legend = c(x = -38, y = 29))
# Did bottom two commented out lines to specify `xyFxd` above  
#xyFxd <- replot()
#dput(xyFxd)


# draw again to get haplotype labels on nodes
x11()
plot(nd3nd4_ntwrk,
  size = sqrt(nd3nd4_sz),
  xy = xyFxd,
  col = "grey70",
  lwd = 2,
  labels = TRUE,
  legend = FALSE)

nd3nd4_hapDF
# Note HaploGroup VII is made up of two FL samples (11 & 12):

##XXX This is the "Atlantic" group separated from rest of Atlantic haplotypes

## Next closest Atlantic population is (sample 13) in Fernandina Beach, FL (group VIII with a VA sample) then in S. Carolina which cluster in group IX, X, and XI

## All rest of FL samples cluster in V (all FL Keys + 1 MS, Texas, & LA = all rest of Parham et al. 2008 Gulf samples except 1) or VI (big bend of FL=the other Parham et al. 2008 Gulf sample)

## 11: Adjacent to BC-30 Spoil Island, Merritt Island, Brevard Co.== just South of St. Augustine/Orlando
## 12: Spoil Island, Merritt Island, Brevard Co.==just South of St. Augustine/Orlando
## 13: Jackson Creek off Amelia River, Nassau Co.==Fernandina Beach area 

################################################################################
# Versions with pasted output
## last updated 2026 July 17
R.Version()[c("platform", "svn rev", "version.string")]
#$platform
#[1] "x86_64-pc-linux-gnu"

#$`svn rev`
#[1] "88974"

#$version.string
#[1] "R version 4.5.2 (2025-10-31)"


packageVersion("ape")
#[1] ‘5.8.1’
packageVersion("pegas")
#[1] ‘1.4’


