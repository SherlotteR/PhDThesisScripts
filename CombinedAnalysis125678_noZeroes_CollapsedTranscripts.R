## Combined Analysis 

## Started 2018/04/12

# About ----------------------------------------------------------------
### This script is to analyse the R18 experiments, replicates 7-10.
### Cosedimentations were performed on 2017.12.11 (7 & 8) and 2017.12.15 (9 & 10)
### All 4 replicates were processed for mass spec from 2017.12.18 - 2017.12.22 and submitted on 2017.12.22
### Mass spec was carried out by Christos Spanos on 2018.04.04 and results received by me on 2018.04.12.
### EDITED ON 2020.04.02 TO COLLAPSE RNA TRANSCRIPTS


# Load data ----------------------------------------------------------------
rm(list = ls())
setwd("\\\\csce.datastore.ed.ac.uk/csce/biology/users/s0906576/Research/Projects/Regulators/14-3-3/R18/LFQ_massSpec")

# Load in csv files where proteins identified with just one peptide have been removed manually (0,1; 1,0; and 1,1)
LFQ1 <- read.csv("E170827/proteinGroups_E170827_Charlotte_LFQ1_no1s.csv")
LFQ2 <- read.csv("E170827/proteinGroups_E170827_Charlotte_LFQ2_no1s.csv")
LFQ5 <- read.csv("E171125/proteinGroups_E171125_Charlotte_LFQ5_no1s.csv")
LFQ6 <- read.csv("E171125/proteinGroups_E171125_Charlotte_LFQ6_no1s.csv")
LFQ7 <- read.csv("E180404/proteinGroups_E180404_Charlotte_7LFQ_no1s.csv")
LFQ8 <- read.csv("E180404/proteinGroups_E180404_Charlotte_8LFQ_updated_no1s.csv")

# Load functions and packages ----------------------------------------------------------------
# 1. Function to remove special characters in regular expressions (needed for grep and sub)
library(stringr)
quotemeta <- function(string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

# 2. plyr 
# for dataframe manipulation, particularly join_all function
library(plyr)

# 3. tidyr
# for data tidying
library(tidyr)

# 4. ggplot2
# for data visualisation
library(ggplot2)

# Data tidying ----------------------------------------------------------------

# Convert Protein.Name columns to character class for searching later
LFQ1$Protein.Name <- as.character(LFQ1$Protein.Name)
LFQ2$Protein.Name <- as.character(LFQ2$Protein.Name)
LFQ5$Protein.Name <- as.character(LFQ5$Protein.Name)
LFQ6$Protein.Name <- as.character(LFQ6$Protein.Name)
LFQ7$Protein.Name <- as.character(LFQ7$Protein.Name)
LFQ8$Protein.Name <- as.character(LFQ8$Protein.Name)



# Pull the simple gene names out of the Protein.Name column (currently they are fasta headers)

  LFQ1$geneName <- lapply(LFQ1$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ1$geneName <- lapply(LFQ1$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ1$geneName <- as.factor(as.character(LFQ1$geneName))

  LFQ2$geneName <- lapply(LFQ2$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ2$geneName <- lapply(LFQ2$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ2$geneName <- as.factor(as.character(LFQ2$geneName))

  LFQ5$geneName <- lapply(LFQ5$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ5$geneName <- lapply(LFQ5$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ5$geneName <- as.factor(as.character(LFQ5$geneName))

  LFQ6$geneName <- lapply(LFQ6$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ6$geneName <- lapply(LFQ6$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ6$geneName <- as.factor(as.character(LFQ6$geneName))

  LFQ7$geneName <- lapply(LFQ7$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ7$geneName <- lapply(LFQ7$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ7$geneName <- as.factor(as.character(LFQ7$geneName))
  
  LFQ8$geneName <- lapply(LFQ8$Protein.Name, gsub, pattern="^.*gn=", replacement = "", ignore.case = TRUE) 
  LFQ8$geneName <- lapply(LFQ8$geneName, gsub, pattern=" .*$", replacement = "", ignore.case = TRUE)
  # if no gene name entered, accession remains
  # and transform that column from a list to a factor because that will just cause problems later
  LFQ8$geneName <- as.factor(as.character(LFQ8$geneName))
  

## Add a 'replicate' column and rename intensity columns
  LFQ1$replicate <- 1
  LFQ1 <- LFQ1[,-2]
  names(LFQ1)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")
  
  LFQ2$replicate <- 2
  names(LFQ2)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")
  
  LFQ5$replicate <- 5
  names(LFQ5)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")
  
  LFQ6$replicate <- 6
  names(LFQ6)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")
  
  LFQ7$replicate <- 7
  names(LFQ7)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")
  
  LFQ8$replicate <- 8
  names(LFQ8)[6:7] <- c("LFQ.intensity.CON", "LFQ.intensity.SAM")



## join the 6 lists together using rbind
LFQall <- rbind(LFQ1[,c(1,4:10)], LFQ2[,c(1,4:10)], LFQ5[,c(1,4:10)], LFQ6[,c(1,4:10)], LFQ7[,c(1,4:10)], LFQ8[,c(1,4:10)])

### duplicate handling 

## collapse RNA transcripts

tran <- grep("-R", LFQall$geneName)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RA$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RB$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RC$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RD$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RE$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- lapply(LFQall$geneName, gsub, pattern="-RG$", replacement = "", ignore.case = TRUE)
LFQall$geneName <- as.factor(as.character(LFQall$geneName))


## check for genes with more than one accession number
length(unique(LFQall$geneName)) #2720
dups <- count(LFQall, vars = "geneName") #gives frequency of each geneName
dups <- subset(dups, dups$freq >= 7) #45


#14-3-3 zeta, Adh, BetaCOP, BetaTub60D, bsf, Cnot4, eEF1delta, Jabba, Lon, Lsd-2, Map205, mei-38, msps, ncd, pyd, sec31, sktl, TER94, tral, zip -- duplicate assignment in mass spec -- remove accession numbers with lower score
LFQall <- LFQall[-c( 1601, #14-3-3 epsilon
                     1526, 3948, 6926, #14-3-3 zeta
                     1574, 3607, 3933, 7127, 8516, 10100, #Adh
                     7216, 10350, #AGO2
                     5969, #akap200
                     3590,  #anxB10
                     1837, 3810, 5571, 7354, 8731, 10431, #betaCOP
                     1785, 3781, 3958, 7363, 8769, 10359,  #betaTub60D
                     1889, 3911, 7452, 8834, 10470, #bsf
                     10424, #CG11148
                     841, 3002, 3942, 6348,  #CG11876
                     3414, 8375, 10133, #CG15618
                     1268, #CG18190
                     1680, 3482, 7292, #CG30122
                     6716, #CG3760
                     7217, #CG4747
                     1902, #CLIP-190
                     1649, 3490, 7436, 8357, 9869, #Cnot4
                     165, 2070, 3547, #coro
                     7093, 8650, 10261, #eEF1delta
                     3934, 3957, 7380, 7381, #Fmr1
                     1556, 3632, #Got1
                     2103, 4154, #GstO2
                     6791, #hts
                     3303, 7192, #Imp
                     8500, 9918, #Jabba
                     3556, 8630, 10212, #Lon
                     1531, 3572, 7095, 10020, #Lsd-2
                     3877, 8830, #Map205
                     1905, 3926, 3927, 3961, 7445, 8855, 10488, #msps
                     1821, 3825, 3951, 7218, 8784, 10383, #ncd
                     6794, #Nopp140
                     836, 3949, 7157, 9274, #Nup98-96
                     3954, 7059, #pcm
                     1133, 1213, 3093, 8565, 8566, 9705, 9790, #pyd
                     10297, #sec31
                     6715, #shep
                     1888, 1907, 3854, 3923, 10329, #shot
                     1644, 3542, 6717, 8108, 9480, #sofe
                     7190, #sqd
                     3888, 10461, #TER94
                     1810, 3959, 7375, 10256,  #tral
                     1877, 3889, 3962, 7433, #Uba1
                     1882, 3908, 3963, 7456, 8723, 10456, #zip
                     3416, 8464, 10082 #zW10
                     ),]


rm(dups)

## Manually add names that are missing
LFQall$geneName <- as.character(LFQall$geneName)

x <- grep("A9UNF9", LFQall$Accession.No)
LFQall[x, 7] <- "Mcm10"

x <- grep("B5RIH0", LFQall$Accession.No)
LFQall[x, 7] <- "CG5815"

x <- grep("K7X538", LFQall$Accession.No)
LFQall[x, 7] <- "armi"

x <- grep("Q6AWL6", LFQall$Accession.No)
LFQall[x, 7] <- "CG8920"

x <- grep("C7LA94", LFQall$Accession.No)
LFQall[x, 7] <- "FI05241p"

x <- grep("Q6IDG3", LFQall$Accession.No)
LFQall[x, 7] <- "Sec10"

x <- grep("K7XI02", LFQall$Accession.No)
LFQall[x, 7] <- "spn-E"

x <- grep("Q4QQC4", LFQall$Accession.No)
LFQall[x, 7] <- "CG2051"

x <- grep("Q7M3J9", LFQall$Accession.No)
LFQall[x, 7] <- "Tubulin beta-3 chain"

x <- grep("K7WSB4", LFQall$Accession.No)
LFQall[x, 7] <- "alphaTub84B"

x <- grep("F6M9W1", LFQall$Accession.No)
LFQall[x, 7] <- "ctrip"

x <- grep("Q9TWZ1", LFQall$Accession.No)
LFQall[x, 7] <- "ERp60"

x <- grep("C7LA94", LFQall$Accession.No)
LFQall[x, 7] <- "Hsp60"

x <- grep("K7WSB4", LFQall$Accession.No)
LFQall[x, 7] <- "alphaTub84B"

x <- grep("Q6NR91", LFQall$Accession.No)
LFQall[x, 7] <- "bt"

x <- grep("Q8T9C4", LFQall$Accession.No)
LFQall[x, 7] <- "SD07683p"

x <- grep("Q6NKL9", LFQall$Accession.No)
LFQall[x, 7] <- "CG1516"

x <- grep("Q4QQC4", LFQall$Accession.No)
LFQall[x, 7] <- "CG2051"

x <- grep("Q95S87", LFQall$Accession.No)
LFQall[x, 7] <- "mRpL53"

x <- grep("Q6AWL6", LFQall$Accession.No)
LFQall[x, 7] <- "CG8920"

x <- grep("Q6IDG3", LFQall$Accession.No)
LFQall[x, 7] <- "Sec10"

x <- grep("Q7M3J7", LFQall$Accession.No)
LFQall[x, 7] <- "msn"

x <- grep("Q7M3J9", LFQall$Accession.No)
LFQall[x, 7] <- "Tubulin beta-3 chain"

x <- grep("F6M9W1", LFQall$Accession.No)
LFQall[x, 7] <- "ctrip"

x <- grep("Q9TWZ1", LFQall$Accession.No)
LFQall[x, 7] <- "ERp60"

x <- grep("E5KZ94", LFQall$Accession.No)
LFQall[x, 7] <- "AGO2"

x <- grep("C7LA94", LFQall$Accession.No)
LFQall[x, 7] <- "Hsp60"

x <- grep("F6M9W1", LFQall$Accession.No)
LFQall[x, 7] <- "ctrip"

x <- grep("Q494I6", LFQall$Accession.No)
LFQall[x, 7] <- "CG1737"

x <- grep("Q4QQC4", LFQall$Accession.No)
LFQall[x, 7] <- "CG2051"

x <- grep("Q7M3J9", LFQall$Accession.No)
LFQall[x, 7] <- "Tubulin beta-3 chain"

x <- grep("Q9TWZ1", LFQall$Accession.No)
LFQall[x, 7] <- "ERp60"

x <- grep("C7LA94", LFQall$Accession.No)
LFQall[x, 7] <- "Hsp60"

x <- grep("K7WSB4", LFQall$Accession.No)
LFQall[x, 7] <- "alphaTub84B"

x <- grep("Q6AWR8", LFQall$Accession.No)
LFQall[x, 7] <- "Psa"

x <- grep("Q5U0W4", LFQall$Accession.No)
LFQall[x, 7] <- "RpL14"

x <- grep("Q4QQC4", LFQall$Accession.No)
LFQall[x, 7] <- "CG2051"

x <- grep("Q6BCZ0", LFQall$Accession.No)
LFQall[x, 7] <- "DppIII"

x <- grep("Q9TWZ1", LFQall$Accession.No)
LFQall[x, 7] <- "ERp60"

x <- grep("Q7M3J9", LFQall$Accession.No)
LFQall[x, 7] <- "Tubulin beta-3 chain"

x <- grep("C7LA94", LFQall$Accession.No)
LFQall[x, 7] <- "Hsp60"

x <- grep("K7WSB4", LFQall$Accession.No)
LFQall[x, 7] <- "alphaTub84B"

x <- grep("A0A0C4FEI6", LFQall$Accession.No)
LFQall[x, 7] <- "Actin-related protein 2/3 complex subunit 3"

x <- grep("A0A0U3DWR7", LFQall$Accession.No)
LFQall[x, 7] <- "Patronin"

x <- grep("A4IJ71", LFQall$Accession.No)
LFQall[x, 7] <- "EndoB"

x <- grep("E4NKI4", LFQall$Accession.No)
LFQall[x, 7] <- "CG18659"

x <- grep("K7WQ31", LFQall$Accession.No)
LFQall[x, 7] <- "armi"

x <- grep("Q29QH4", LFQall$Accession.No)
LFQall[x, 7] <- "CG14036"

x <- grep("Q29QV6", LFQall$Accession.No)
LFQall[x, 7] <- "FANCI"

x <- grep("Q6GKZ3", LFQall$Accession.No)
LFQall[x, 7] <- "HDAC4"

x <- grep("Q95SD0", LFQall$Accession.No)
LFQall[x, 7] <- "GM03203p"


## Remove the columns we don't need 
LFQnew <- LFQall[,c("geneName", "LFQ.intensity.SAM", "LFQ.intensity.CON", 
                    "replicate", "Accession.No")]



# Calculations -------------------------------------------------------------------
## Make zeroes NA
  
LFQnew$LFQ.intensity.SAM <- as.numeric(sub("^0$", NA, LFQnew$LFQ.intensity.SAM)) #substitute zeros for NA
LFQnew$LFQ.intensity.CON <- as.numeric(sub("^0$", NA, LFQnew$LFQ.intensity.CON)) #substitute zeros for NA
  
  
## Next, calculate log(2) for each intensity
  
LFQnew$logSAM <- log2(LFQnew$LFQ.intensity.SAM)
LFQnew$logCON <- log2(LFQnew$LFQ.intensity.CON)



## t-tests -------------------------------------------------------------------
### paired, un-normalised ----------------------------------------------------------------------
## Carry out paired t-test for each protein, using the log.intensity value

onePair <- NULL
for (i in unique(LFQnew$geneName)) {
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, LFQall$geneName)
  a <- LFQnew[rows,]
  comp <- complete.cases(a[,c(2:3)])
  l <- length(subset(comp, comp == TRUE))
  if (l >= 2){
    tt <- with(a, t.test(logSAM, logCON, paired = TRUE))
    LFQnew[rows, "pairedp.value"] <- tt$p.value
  } else {
    onePair[rows] <- LFQnew[rows, "geneName"]
    LFQnew[rows, "pairedp.value"] <- NA
  }
  
}

onePair <- na.omit(onePair)
write.csv(onePair, "onePairR18Names.csv")

## Calculate ratios. Proteins that don't have at least 1 SAM value and 1 CON value are not in this dataset!

d <- NULL
for (i in unique(LFQnew$geneName)){
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, LFQnew$geneName)
  b <- LFQnew[rows,]
  a <-  b[complete.cases(b[,c(2:3)]) == TRUE,]
  if (nrow(a) >= 1){
    if ((nrow(a[a$replicate == 1,])) >= 1){
      a$logSAM1 <- a[a$replicate == 1, 6]
      a$logCON1 <- a[a$replicate == 1, 7]
      
    } else {
      a$logSAM1 <- NA
      a$logCON1 <- NA
    }
    if ((nrow(a[a$replicate == 2,])) >= 1){
      a$logSAM2 <- a[a$replicate == 2, 6]
      a$logCON2 <- a[a$replicate == 2, 7]
      
    } else {
      a$logSAM2 <- NA
      a$logCON2 <- NA
    }
    if ((nrow(a[a$replicate == 5,])) >= 1){
      a$logSAM5 <- a[a$replicate == 5, 6]
      a$logCON5 <- a[a$replicate == 5, 7]
      
    } else {
      a$logSAM5 <- NA
      a$logCON5 <- NA
    }
    if ((nrow(a[a$replicate == 6,])) >= 1){
      a$logSAM6 <- a[a$replicate == 6, 6]
      a$logCON6 <- a[a$replicate == 6, 7]
      
    } else {
      a$logSAM6 <- NA
      a$logCON6 <- NA
    }
    
    if ((nrow(a[a$replicate == 7,])) >= 1){
      a$logSAM7 <- a[a$replicate == 7, 6]
      a$logCON7 <- a[a$replicate == 7, 7]
      
    } else {
      a$logSAM7 <- NA
      a$logCON7 <- NA
    }
    if ((nrow(a[a$replicate == 8,])) >= 1){
      a$logSAM8 <- a[a$replicate == 8, 6]
      a$logCON8 <- a[a$replicate == 8, 7]
      
    } else {
      a$logSAM8 <- NA
      a$logCON8 <- NA
      
    }
    
    a$logSAM <- mean(a$logSAM, na.rm = TRUE)
    a$logCON <- mean(a$logCON, na.rm = TRUE)
    a$logratio <- a$logSAM - a$logCON
    
    d <- rbind(d, a[1,])
  }
}

setwd("\\\\csce.datastore.ed.ac.uk/csce/biology/users/s0906576/Research/Projects/Regulators/14-3-3/R18/LFQ_massSpec/R18_125678_complexesAndInteractionsAnalysis/R18Script_outputs")
write.csv(d, "d_fromCombAnal125678_noZero_collRNA.csv")     
xx <- subset(d, d$logratio > 1 & d$pairedp.value < 0.01)

library(see)
## this is how I plotted for my poster but this is different data
p1<-ggplot(data = d, aes(x = logratio, y = -log10(pairedp.value), text=geneName))+
  theme_light(base_size = 15, base_family = "")+
  geom_point(colour = "grey45", size = 1.5, alpha = 0.8, shape = 16) +
  scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "Fold Change in Microtubule Binding", y = "p value")+
# geom_point(data = subset(d, d$logSAM > 23))+
# geom_text(data = xx, label = xx$geneName, nudge_x = 0.17, size = 3)
  geom_point(data = subset(d, d$geneName == "ncd"), colour = "red", size = 2)+
  geom_point(data = subset(d, d$geneName == "pav"), colour = "yellow2", size = 2)+
  geom_point(data = subset(d, d$geneName == "tum"), colour = "orange", size = 2)+
  geom_point(data = subset(d, d$geneName == "Incenp"), colour = "blue", size = 2)+
  #geom_point(data = subset(d, d$geneName == "borr"), colour = "blue", size = 2)+
  geom_point(data = subset(d, d$geneName == "ial"), colour = "blue", size = 2)#+
  #geom_point(data = subset(d, d$geneName == "sle"), colour = "green2", size = 2)+
  #geom_point(data = subset(d, d$geneName == "stai"), colour = "green2", size = 2)+
  #geom_point(data = subset(d, d$geneName == "Rm62"), colour = "green2", size = 2)+
  #geom_point(data = subset(d, d$geneName == "CG12909"), colour = "green2", size = 2)+
  #geom_point(data = subset(d, d$geneName == "Ns3"), colour = "green2", size = 2)+
  #geom_point(data = subset(d, d$geneName == "ast"), colour = "blue", size = 2)+
  #geom_point(data = subset(d, d$geneName == "l(3)07882"), colour = "green2", size = 2)

library(plotly)
p<-ggplotly(p1, tooltip = 'text')
Sys.setenv("plotly_username"="#########")
Sys.setenv("plotly_api_key"="#####################")
api_create(p, filename = "R18_interactiveGraph")
  
  
ggplot(data = d, aes(x = logratio, y = -log10(pairedp.value)))+
  theme_light(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "", y = "p value")+
  geom_point(data = subset(d, d$logSAM > 23))+
  #geom_label(data = z, label = geneName)+
  geom_point(data = subset(d, d$geneName == "Bruce"), colour = "red", size = 2)
  #geom_point(data = subset(d, d$logratio > "1"), colour = "red", size = 2)

spliceosome <- z[z$geneName %in% c("CG10333", "Cdc5", "Prp8", "Prp19", "Prp31", "Rm62", "SmD2", "B52"),]
MTcytoOrg <- z[z$geneName %in% c("CG11120", "ncd", "tum", "borr", "pav", "Incenp", "ial", "mip130", "Chro", "aPKC", "SmD2", "Fib", "larp", "RpS16", "Ns1"),]
mitSpinAss <- z[z$geneName %in% c("ncd", "tum", "borr", "pav", "Incenp", "ial", "mip130", "Ns1"),]

ggplot(data = d, aes(x = logratio, y = -log10(pairedp.value)))+
  theme_gray(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "", y = "p value")+
  geom_point(data = subset(d, d$logSAM > 23))+
  #geom_label(data = z, label = geneName)+
  geom_point(data = spliceosome, colour = "yellow", size = 2)+
  geom_point(data = MTcytoOrg, colour = "red", size = 2)+
  geom_point(data = mitSpinAss, colour = "blue", size = 2)+
  geom_point(data = subset(spliceosome, spliceosome$geneName == "SmD2"), colour = "yellow", size = 2)+
  geom_point(data = subset(spliceosome, spliceosome$geneName == "SmD2"), colour = "red", size = 1)

#write.csv(d[order(d$pairedp.value),], "XXX") #require 4 or more con/sam values for inclusion in t-test
#write.csv(d[order(d$pairedp.value),], "XXX") #require 6 or more
#write.csv(d[order(d$pairedp.value),], "XXX") #require 8 or more

LFQnew$rat <- LFQnew$logSAM - LFQnew$logCON

LFQnew$replicate <- as.factor(LFQnew$replicate)
LFQnew$ftt <- NA
levl1433 <- (LFQnew[grep('14-3-3epsilon', LFQnew$geneName),])
levl1433 <- levl1433[-1,]
LFQnew[LFQnew$replicate == "1", 'ftt' ] <- levl1433[1, 'rat']
LFQnew[LFQnew$replicate == "2", 'ftt' ] <- levl1433[2, 'rat']
LFQnew[LFQnew$replicate == "5", 'ftt' ] <- levl1433[3, 'rat']
LFQnew[LFQnew$replicate == "6", 'ftt' ] <- levl1433[4, 'rat']
LFQnew[LFQnew$replicate == "7", 'ftt' ] <- levl1433[5, 'rat']
LFQnew[LFQnew$replicate == "8", 'ftt' ] <- levl1433[6, 'rat']

ggplot(data=LFQnew, aes(x = LFQnew$ftt, y = LFQnew$rat))+
  geom_point()


install.packages("broom")
library(broom)

mod2 <- glm(LFQnew$rat ~ LFQnew$ftt/LFQnew$geneName)
options(max.print=1000000)
hits <- summary(mod2)
tidy_mod2 <- tidy(mod2)
write.csv(tidy_mod2, "hits.csv")
plot(mod2)


## no effect of 14-3-3epsilon levels on overall protein content on microtubules (despite significance in model, there are not 7018 dfs in this model so I don't believe the stars)

kwn.MAPs <- c("ncd", "Eb1", "tacc", "dgt2", "dgt3", "dgt4", "dgt5", "dgt6", "msps", "ssp2", 
              "sub", "pav", "KLP61F", "KLP54D", "ASP", "Nod", "ial", "KLP10A", "Incenp", 
              "Mars", "Mei-38", "Grip71", "Cyclin B", "Mud", "Cnn", "Feo", "Axs", "chb", 
              "Polo", "KLP67A", "KLP3A")
e <- c(rep(NA, 10))
names(e) <- names(LFQnew)
for (i in kwn.MAPs){
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rowNo <- grep(searchterm, LFQnew$geneName, ignore.case = TRUE)
  e <- rbind(e, LFQnew[rowNo,])
}

e <- e[-1, ]

mod3 <- glm(e$rat ~ e$ftt/e$geneName)
summary(mod3)
plot(mod3)


ggplot(data=e, aes(x = e$ftt, y = e$rat))+
  geom_point()+
  geom_point(aes(x = e[1, 'ftt'], y = e[1, 'rat']), colour = "red")+
  geom_point(aes(x = e[2, 'ftt'], y = e[2, 'rat']), colour = "red")+
  geom_point(aes(x = e[3, 'ftt'], y = e[3, 'rat']), colour = "red")+
  geom_point(aes(x = e[4, 'ftt'], y = e[4, 'rat']), colour = "red")+
  geom_point(aes(x = e[5, 'ftt'], y = e[5, 'rat']), colour = "red")+
  geom_point(aes(x = e[6, 'ftt'], y = e[6, 'rat']), colour = "red")+
  geom_abline(slope = 0.673870189, intercept = 0.125957588)

ggplot(data=e, aes(x = e$ftt, y = e$rat))+
  geom_point()+
  geom_point(aes(x = e[51, 'ftt'], y = e[51, 'rat']), colour = "red")+
  geom_point(aes(x = e[52, 'ftt'], y = e[52, 'rat']), colour = "red")+
  geom_point(aes(x = e[53, 'ftt'], y = e[53, 'rat']), colour = "red")+
  geom_point(aes(x = e[54, 'ftt'], y = e[54, 'rat']), colour = "red")+
  geom_point(aes(x = e[55, 'ftt'], y = e[55, 'rat']), colour = "red")+
  geom_point(aes(x = e[56, 'ftt'], y = e[56, 'rat']), colour = "red")

write.csv(e, "knownMAPSratios2018.0416.csv")
### ---------------------------------------

#correct ratios by 14-3-3 levels
levl1433 <- (LFQnew[grep('14-3-3epsilon', LFQnew$geneName),])
levl1433 <- levl1433[-1,]
LFQnew[LFQnew$replicate == "1", 'fttLogSAM' ] <- levl1433[1, 'LFQ.intensity.SAM']
LFQnew[LFQnew$replicate == "2", 'fttLogSAM' ] <- levl1433[2, 'LFQ.intensity.SAM']
LFQnew[LFQnew$replicate == "5", 'fttLogSAM' ] <- levl1433[3, 'LFQ.intensity.SAM']
LFQnew[LFQnew$replicate == "6", 'fttLogSAM' ] <- levl1433[4, 'LFQ.intensity.SAM']
LFQnew[LFQnew$replicate == "7", 'fttLogSAM' ] <- levl1433[5, 'LFQ.intensity.SAM']
LFQnew[LFQnew$replicate == "8", 'fttLogSAM' ] <- levl1433[6, 'LFQ.intensity.SAM']

LFQnew[LFQnew$replicate == "1", 'fttLogCON' ] <- levl1433[1, 'LFQ.intensity.CON']
LFQnew[LFQnew$replicate == "2", 'fttLogCON' ] <- levl1433[2, 'LFQ.intensity.CON']
LFQnew[LFQnew$replicate == "5", 'fttLogCON' ] <- levl1433[3, 'LFQ.intensity.CON']
LFQnew[LFQnew$replicate == "6", 'fttLogCON' ] <- levl1433[4, 'LFQ.intensity.CON']
LFQnew[LFQnew$replicate == "7", 'fttLogCON' ] <- levl1433[5, 'LFQ.intensity.CON']
LFQnew[LFQnew$replicate == "8", 'fttLogCON' ] <- levl1433[6, 'LFQ.intensity.CON']


LFQnew$corrLogSAM <- LFQnew$LFQ.intensity.SAM * LFQnew$fttLogSAM
LFQnew$corrLogCON <- LFQnew$LFQ.intensity.CON * LFQnew$fttLogCON


# t-tests
# paired, un-normalised 
## Carry out paired t-test for each protein, using the log.intensity value

for (i in unique(LFQnew$geneName)) {
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, LFQall$geneName)
  a <- LFQnew[rows,]
  comp <- complete.cases(a[,c(2:3)])
  l <- length(subset(comp, comp == TRUE))
  if (l >= 2){
    tt <- with(a, t.test(corrLogSAM, corrLogCON, paired = TRUE))
    LFQnew[rows, "pairedp.value"] <- tt$p.value
  } else {
    LFQnew[rows, "pairedp.value"] <- NA
  }
  
}

LFQnew$corrLogSAM <- log(LFQnew$corrLogSAM, 2)
LFQnew$corrLogCON <- log(LFQnew$corrLogCON, 2)


## Calculate ratios

f <- NULL
for (i in unique(LFQnew$geneName)){
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, LFQnew$geneName)
  b <- LFQnew[rows,]
  a <-  b[complete.cases(b[,c(2:3)]) == TRUE,]
  if (nrow(a) >= 1){
    if ((nrow(a[a$replicate == 1,])) >= 1){
      a$logcorrSAM1 <- a[a$replicate == 1, 12]
      a$logcorrCON1 <- a[a$replicate == 1, 13]
    } else {
      a$logcorrSAM1 <- NA
      a$logcorrCON1 <- NA
    }
    if ((nrow(a[a$replicate == 2,])) >= 1){
      a$logcorrSAM2 <- a[a$replicate == 2, 12]
      a$logcorrCON2 <- a[a$replicate == 2, 13]
      
    } else {
      a$logcorrSAM2 <- NA
      a$logcorrCON2 <- NA
    }
    if ((nrow(a[a$replicate == 5,])) >= 1){
      a$logcorrSAM5 <- a[a$replicate == 5, 12]
      a$logcorrCON5 <- a[a$replicate == 5, 13]
      
    } else {
      a$logcorrSAM5 <- NA
      a$logcorrCON5 <- NA
    }
    if ((nrow(a[a$replicate == 6,])) >= 1){
      a$logcorrSAM6 <- a[a$replicate == 6, 12]
      a$logcorrCON6 <- a[a$replicate == 6, 13]
      
    } else {
      a$logcorrSAM6 <- NA
      a$logcorrCON6 <- NA
    }
    
    if ((nrow(a[a$replicate == 7,])) >= 1){
      a$logcorrSAM7 <- a[a$replicate == 7, 12]
      a$logcorrCON7 <- a[a$replicate == 7, 13]
      
    } else {
      a$logcorrSAM7 <- NA
      a$logcorrCON7 <- NA
    }
    if ((nrow(a[a$replicate == 8,])) >= 1){
      a$logcorrSAM8 <- a[a$replicate == 8, 12]
      a$logcorrCON8 <- a[a$replicate == 8, 13]
      
    } else {
      a$logcorrSAM8 <- NA
      a$logcorrCON8 <- NA
      
    }
    
    a$logcorrSAM <- mean(a$corrLogSAM, na.rm = TRUE)
    a$logcorrCON <- mean(a$corrLogCON, na.rm = TRUE)
    a$logcorrratio <- a$logcorrSAM - a$logcorrCON
    
    f <- rbind(f, a[1,])
  }
}


ggplot(data = f, aes(x = logcorrratio, y = -log10(pairedp.value)))+
  theme_light(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "", y = "p value")+
  geom_point(data = subset(f, f$logcorrSAM > 23))+
  geom_point(data = subset(f, f$logcorrratio > "1"), colour = "red", size = 2)


## --------------------------------------------
# Deal with those proteins present in SAM but not CON
## put fake zero numbers into CON column (z.2)

y <- LFQnew[is.na(LFQnew$pairedp.value),]
y <- y[is.na(y$logCON),]
y <- y[is.na(y$logSAM) == FALSE, ]

z <- data.frame()
newrow <- NULL
for (i in unique(y$geneName)){
  searchterm1 <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm1, y$geneName)
  b <- LFQnew[rows,]
  geneName <- b$geneName
  if (nrow(b)>1){
    logSAM <- mean(b$logSAM)
  } else {
    logSAM <- b$logSAM
  }
  sdSAM <- sd(b$logSAM)
  no.reps <- nrow(b)
  newrow <- data.frame(geneName, no.reps, logSAM, sdSAM)
  z <- rbind(z, newrow)
}

plot(z$logSAM)

z.1 <- subset(z, z$logSAM >23)

z.2 <- read.csv("NoCONproteins20180503.csv")

for (i in unique(z.2$geneName)) {
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, z.2$geneName)
  a <- z.2[rows,]
  comp <- complete.cases(a[,c(2:3)])
  l <- length(subset(comp, comp == TRUE))
  if (l >= 2){
    tt <- with(a, t.test(logSAM, logCON, paired = TRUE))
    z.2[rows, "pairedp.value"] <- tt$p.value
  } else {
    z.2[rows, "pairedp.value"] <- NA
  }
  
}

z.d <- NULL
for (i in unique(z.2$geneName)){
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, z.2$geneName)
  a <- z.2[rows,]
  if (nrow(a) >= 1){
    if ((nrow(a[a$replicate == 1,])) >= 1){
      a$logSAM1 <- a[a$replicate == 1, 7]
      a$logCON1 <- a[a$replicate == 1, 8]
      
    } else {
      a$logSAM1 <- NA
      a$logCON1 <- NA
    }
    if ((nrow(a[a$replicate == 2,])) >= 1){
      a$logSAM2 <- a[a$replicate == 2, 7]
      a$logCON2 <- a[a$replicate == 2, 8]
      
    } else {
      a$logSAM2 <- NA
      a$logCON2 <- NA
    }
    if ((nrow(a[a$replicate == 5,])) >= 1){
      a$logSAM5 <- a[a$replicate == 5, 7]
      a$logCON5 <- a[a$replicate == 5, 8]
      
    } else {
      a$logSAM5 <- NA
      a$logCON5 <- NA
    }
    if ((nrow(a[a$replicate == 6,])) >= 1){
      a$logSAM6 <- a[a$replicate == 6, 7]
      a$logCON6 <- a[a$replicate == 6, 8]
      
    } else {
      a$logSAM6 <- NA
      a$logCON6 <- NA
    }
    
    if ((nrow(a[a$replicate == 7,])) >= 1){
      a$logSAM7 <- a[a$replicate == 7, 7]
      a$logCON7 <- a[a$replicate == 7, 8]
      
    } else {
      a$logSAM7 <- NA
      a$logCON7 <- NA
    }
    if ((nrow(a[a$replicate == 8,])) >= 1){
      a$logSAM8 <- a[a$replicate == 8, 7]
      a$logCON8 <- a[a$replicate == 8, 8]
      
    } else {
      a$logSAM8 <- NA
      a$logCON8 <- NA
      
    }
    
    a$logSAM <- mean(a$logSAM, na.rm = TRUE)
    a$logCON <- mean(a$logCON, na.rm = TRUE)
    a$logratio <- a$logSAM - a$logCON
    
    z.d <- rbind(z.d, a[1,])
  }
}

## This little loop creates the NoCON20201024.csv list. It checks which proteins from Z (proteins with at least 1 NA in CON columns) are actually present in d (e.g. plotted, and therefore not NoCON) and returns just the part of z that is unique

zo <-list()
zoe <- as.numeric(NA)
for (i in z$geneName){ 
  zoe<-grep(i, d$geneName)
  if(length(zoe)== 0){
    zo <- c(zo, i)
  }
}
#write.csv(subset(z, z$geneName %in% zo), "NoCONproteins_20201024.csv")

ggplot(data = z.d, aes(x = logratio, y = -log10(pairedp.value)))+
  theme_light(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  #scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  #scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "", y = "p value")+
  geom_point(data = subset(z.d, z.d$geneName == "M1BP"), colour = "red", size = 2)+
  geom_point(data = subset(z.d, z.d$geneName == "Det"), colour = "skyblue2", size = 2)
    
#--------------------------------------------------
#comparing to robin's pulldowns

bindingPros <- read.csv("RobinPulldown_1433BindingProteins.csv")
p1 + geom_point(data = subset(d,d$geneName %in% bindingPros$RobinGST == TRUE), colour = "yellow")
p1 + geom_point(data = subset(d,d$geneName %in% bindingPros$RobinMBP == TRUE), colour = "yellow")
p1 + geom_point(data = subset(d,d$geneName %in% bindingPros$RobinS2 == TRUE), colour = "yellow")

#known MAPs
kwn.MAPs <- c("ncd", "Eb1", "tacc", "dgt2", "dgt3", "dgt4", "dgt5", "dgt6", "msps", "ssp2", 
              "sub", "pav", "Klp61F", "Klp54D", "asp-RA", "Nod", "ial", "Klp10A", "Incenp", 
              "mars", "mei-38", "Grip71", "Cyclin B", "Mud", "cnn", "feo", "Axs", "chb", 
              "polo", "Klp67A", "Klp3A")
#suspected MAPs
spt.MAPs <- c("Ndc80", "Nuf", "Grip84", "CG6015", "Chro", "Sticky", "CG1951", "SAK kinase (Plk4)",
              "Map205", "awd", "gwl", "cmet", "spindly")

p2<-p1 + geom_point(data = subset(d,d$geneName %in% kwn.MAPs == TRUE), colour = "green4")
p2<-p2 + geom_point(data = subset(d,d$geneName %in% spt.MAPs == TRUE), colour = "green2")
p2

library(plotly)
p<-ggplotly(p2, tooltip = 'text')
Sys.setenv("plotly_username"="C.Repton")
Sys.setenv("plotly_api_key"="Qoixsc6k4EgXg34oG5SY")
api_create(p, filename = "R18_interactiveGraphKwnSptMAPs")

MTBinding <- read.csv("MTBindingProteins_Flybase.csv")
femSpindle <- c("nod", "gammaTub37C", "klp10A", "Cks30A", "grip71")
polar <- read.csv("oocytePolarity.csv")
splice <- read.csv("RNAsplicing.csv")

p3<-p1 + geom_point(data = subset(d,d$geneName %in% MTBinding$SYMBOL == TRUE), colour = "green4")
p3<-p3 + geom_point(data = subset(d,d$geneName %in% femSpindle == TRUE), colour = "blue4")
p3<-p3 + geom_point(data = subset(d,d$geneName %in% polar$SYMBOL == TRUE), colour = "red4")
p3<-p3 + geom_point(data = subset(d,d$geneName %in% splice$SYMBOL == TRUE), colour = "gold4")
p3


ggplot(data = d, aes(x = logratio, y = -log10(pairedp.value)))+
  theme_light(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  geom_point()+
  gghighlight(logratio > 0.5, pairedp.value < 0.05)+
  scale_x_continuous(limits = c(-2, 3), labels = c("1/4", "1/2", "1x", "2x", "4x", "8x"))+
  scale_y_continuous(limits = c(0, 4), labels = c("1", "0.1", "0.01", "0.001", "0.0001"))+
  labs(x = "", y = "p value")
#--------------------------------------------
#Make file of one pair or fewer proteins
OP <- NULL
for (i in unique(onePair)){
  searchterm <- paste("^", quotemeta(i), "$", sep = "")
  rows <- grep(searchterm, LFQnew$geneName)
  a <- data.frame(matrix(NA, nrow = 6, ncol = 8))
  names(a) <- c(names(LFQnew))
  a[1,] <- LFQnew[rows[1],]
  a[2,] <- LFQnew[rows[2],]
  a[3,] <- LFQnew[rows[3],]
  a[4,] <- LFQnew[rows[4],]
  a[5,] <- LFQnew[rows[5],]
  a[6,] <- LFQnew[rows[6],]
  
  if(length(na.omit(a[a$replicate == 1, 6])) >= 1){
    a$logSAM1[1] <- a[a$replicate == 1, 6][1]
  } else {
    a$logSAM1[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 2, 6])) >= 1){
    a$logSAM2[1] <- a[a$replicate == 2, 6][1]
  } else {
    a$logSAM2[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 5, 6])) >= 1){
    a$logSAM5[1] <- a[a$replicate == 5, 6][1]
  } else {
    a$logSAM5[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 6, 6])) >= 1){
    a$logSAM6[1] <- a[a$replicate == 6, 6][1]
  } else {
    a$logSAM6[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 7, 6])) >= 1){
    a$logSAM7[1] <- a[a$replicate == 7, 6][1]
  } else {
    a$logSAM7[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 8, 6])) >= 1){
    a$logSAM8[1] <- a[a$replicate == 8, 6][1]
  } else {
    a$logSAM8[1] <- NA
  }
  
  
  if(length(na.omit(a[a$replicate == 1, 7])) >= 1){
    a$logCON1[1] <- a[a$replicate == 1, 7][1]
  } else {
    a$logCON1[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 2, 7])) >= 1){
    a$logCON2[1] <- a[a$replicate == 2, 7][1]
  } else {
    a$logCON2[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 5, 7])) >= 1){
    a$logCON5[1] <- a[a$replicate == 5, 7][1]
  } else {
    a$logCON5[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 6, 7])) >= 1){
    a$logCON6[1] <- a[a$replicate == 6, 7][1]
  } else {
    a$logCON6[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 7, 7])) >= 1){
    a$logCON7[1] <- a[a$replicate == 7, 7][1]
  } else {
    a$logCON7[1] <- NA
  }
  if(length(na.omit(a[a$replicate == 8, 7])) >= 1){
    a$logCON8[1] <- a[a$replicate == 8, 7][1]
  } else {
    a$logCON8[1] <- NA
  }
  
  a$logSAM <- mean(a$logSAM, na.rm = TRUE)
  a$logCON <- mean(a$logCON, na.rm = TRUE)
  a$logratio <- a$logSAM - a$logCON
  a$numSAMpos <- sum(!is.na(a[1, 9:14])) 
  a$numCONpos <- sum(!is.na(a[1, 15:20])) 
  a$numPosTotal <- a$numCONpos + a$numSAMpos
  OP <- rbind(OP, a[1,])
}

write.csv(OP, "onePairProteins.csv")

##minimal version of graph for 2nd year poster (no labels, custom colours, etc)
p10<-ggplot(data = d, aes(x = logratio, y = -log10(pairedp.value), text=geneName))+
  theme_classic(base_size = 15, base_family = "")+
  geom_point(colour = "grey") +
  scale_x_continuous(limits = c(-2, 3), labels = c("", "", "", "", "", ""))+
  scale_y_continuous(limits = c(0, 4), labels = c("", "", "", "", ""))+
  labs(x = "", y = "")+
  # geom_text(data = xx, label = xx$geneName, nudge_x = 0.17, size = 3)+
  annotate("rect", xmin=0.72, xmax=1.7, ymin=1.47, ymax=4, alpha = "0.2", fill = "#a2cc82")+
  geom_point(data = subset(d, d$geneName == "ncd"), colour = "#be1622", size = 3)+
  geom_point(data = subset(d, d$geneName == "pav"), colour = "#f39200", size = 3)+
  geom_point(data = subset(d, d$geneName == "tum"), colour = "#f39200", size = 3)+
  geom_point(data = subset(d, d$geneName == "Incenp"), colour = "#1d71b8", size = 3)+
  geom_point(data = subset(d, d$geneName == "borr"), colour = "#1d71b8", size = 3)+
  geom_point(data = subset(d, d$geneName == "ial"), colour = "#1d71b8", size = 3)+
  geom_point(data = subset(d, d$geneName == "sle"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "stai"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "Rm62"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "CG12909"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "Ns3"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "ast"), colour = "#a2cc82", size = 3)+
  geom_point(data = subset(d, d$geneName == "l(3)07882"), colour = "#a2cc82", size = 3)

## --------------------------------------------------------------------------------------
# 20200318 Looking for complexes in the data
## Are there interactors binding microtubules in this data? Do they change together, in the presence of R18?

## I want to look at positive and negative sides of the graph, but I'll keep them separate. Previously I have set the ratio cutoff
## to 1.5x but I think I'll use 1.1x (= 0.15 in log2) today.

## subset the dataset d into significant (below 0.1 in this case) points, and then positive/negative
Significant <- subset(d, d$pairedp.value < 0.1)
IncreaseMT <- subset(Significant, Significant$logratio > 0.15)
DecreaseMT <- subset(Significant, Significant$logratio < -0.15)

## export IncreaseMT and DecreaseMT to look at in databases online
write.csv(IncreaseMT, "./R18_125678_complexesAndInteractionsAnalysis/R18Script_outputs/R18_0.1pvalue_0.15logratio.csv")
write.csv(DecreaseMT, "./R18_125678_complexesAndInteractionsAnalysis/R18Script_outputs/R18_0.1pvalue_-0.15logratio.csv")


## make a subset of those with zeroes -- no p-value?
naPvalue <- subset(d, is.na(d$pairedp.value) == TRUE)
IncreaseMT.naP <- subset(naPvalue, naPvalue$logratio > 0.15)
DecreaseMT.naP <- subset(naPvalue, naPvalue$logratio < -0.15)

write.csv(IncreaseMT.naP, "./R18_125678_complexesAndInteractionsAnalysis/R18Script_outputs/R18_NApvalue_0.15logratio.csv")
write.csv(DecreaseMT.naP, "./R18_125678_complexesAndInteractionsAnalysis/R18Script_outputs/R18_NApvalue_-0.15logratio.csv")

#--------------------------------------------------
#thesis
## make jitter plot to compare no con proteins to sam + con
library(ggplot2)
n <- read.csv("NoCONproteins_20201024.csv")

compd <- d[, c(1, 6)]
compd$p <- "Detected in controls"
compn <- subset(n, n$no.reps > 1)
compn <- subset(compn, is.na(compn$logSAM) == FALSE)
compn <- compn[, c(2,4)]
compn$p <- "Not detected in controls"
comp <- rbind(compd, compn)

median(compd$logSAM) #21.96466
median(compn$logSAM) #19.45807

p2<-ggplot(comp, aes(y = logSAM, x = p))+
  geom_jitter()+
  geom_segment(x= 0.7 , y = 21.96466, xend = 1.3 , yend = 21.96466, colour = "red", size = 1.2)+
  geom_segment(x= 1.7 , y = 19.45807, xend = 2.3 , yend = 19.45807, colour = "red", size = 1.2)+
  labs(x = "", y = "Sample Intensity, millions")+
  scale_y_continuous(breaks = c(19.93157, 23.2535, 26.57542, 29.89735), labels = c("1", "10", "100", "1000"))+
  theme_bw()+
  theme(text = element_text(size = 15))

### Thesis
### Making sam vs con graphs for all replicates (including ones I didn't end up using)

LFQ3 <- read.csv("E170827/proteinGroups_E170827_Charlotte_LFQ0_no1s.csv")
LFQ4 <- read.csv("E171125/proteinGroups_E171125_Charlotte_LFQ4_no1s.csv")
LFQ9 <- read.csv("E180404/proteinGroups_E180404_Charlotte_9LFQ_no1s.csv")
LFQ10 <- read.csv("E180404/proteinGroups_E180404_Charlotte_10LFQ_no1s.csv")

LFQ3$LFQ.intensity.SAM <- as.numeric(sub("^0$", NA, LFQ3$LFQ.intensity.SAM)) #substitute zeros for NA
LFQ3$LFQ.intensity.CON <- as.numeric(sub("^0$", NA, LFQ3$LFQ.intensity.CON)) #substitute zeros for NA
LFQ4$LFQ.intensity.S4 <- as.numeric(sub("^0$", NA, LFQ4$LFQ.intensity.S4)) #substitute zeros for NA
LFQ4$LFQ.intensity.C4 <- as.numeric(sub("^0$", NA, LFQ4$LFQ.intensity.C4)) #substitute zeros for NA
LFQ9$LFQ.intensity.SAM <- as.numeric(sub("^0$", NA, LFQ9$LFQ.intensity.SAM)) #substitute zeros for NA
LFQ9$LFQ.intensity.CON <- as.numeric(sub("^0$", NA, LFQ9$LFQ.intensity.CON)) #substitute zeros for NA
LFQ10$LFQ.intensity.SAM <- as.numeric(sub("^0$", NA, LFQ10$LFQ.intensity.SAM)) #substitute zeros for NA
LFQ10$LFQ.intensity.CON <- as.numeric(sub("^0$", NA, LFQ10$LFQ.intensity.CON)) #substitute zeros for NA

## Next, calculate log(2) for each intensity

LFQ3$logSAM <- log2(LFQ3$LFQ.intensity.SAM)
LFQ3$logCON <- log2(LFQ3$LFQ.intensity.CON)
LFQ4$logS4 <- log2(LFQ4$LFQ.intensity.S4)
LFQ4$logC4 <- log2(LFQ4$LFQ.intensity.C4)
LFQ9$logSAM <- log2(LFQ9$LFQ.intensity.SAM)
LFQ9$logCON <- log2(LFQ9$LFQ.intensity.CON)
LFQ10$logSAM <- log2(LFQ10$LFQ.intensity.SAM)
LFQ10$logCON <- log2(LFQ10$LFQ.intensity.CON)


g1<-ggplot(data = d, aes(logSAM1, logCON1))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 1", y = "Log Control 1")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g2<-ggplot(data = d, aes(logSAM2, logCON2))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 2", y = "Log Control 2")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g3<-ggplot(data = LFQ3, aes(logSAM, logCON))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 3", y = "Log Control 3")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g4<-ggplot(data = LFQ4, aes(logS4, logC4))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 4", y = "Log Control 4")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g5<-ggplot(data = d, aes(logSAM5, logCON5))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 5", y = "Log Control 5")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g6<-ggplot(data = d, aes(logSAM6, logCON6))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 6", y = "Log Control 6")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g7<-ggplot(data = d, aes(logSAM7, logCON7))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 7", y = "Log Control 7")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g8<-ggplot(data = d, aes(logSAM8, logCON8))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 8", y = "Log Control 8")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g9<-ggplot(data = LFQ9, aes(logSAM, logCON))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 9", y = "Log Control 9")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

g10<-ggplot(data = LFQ10, aes(logSAM, logCON))+
  theme_light(base_size = 15, base_family = "")+
  labs(x = "Log Sample 10", y = "Log Control 10")+
  scale_x_continuous(limits = c(15,35))+
  scale_y_continuous(limits = c(15,35))+
  geom_point(size = 1)+
  geom_abline(aes(slope = 1, intercept = 0))

library("cowplot")
plot_grid(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, ncol = 4, nrow = 3)
rsq <- function(x, y) summary(lm(y~x))$r.squared
rsq(d$logSAM1, d$logCON1) #0.871
rsq(d$logSAM2, d$logCON2) #0.957
rsq(LFQ3$logSAM, LFQ3$logCON) #0.520
rsq(LFQ4$logS4, LFQ4$logC4) #0.370
rsq(d$logSAM5, d$logCON5) #0.852
rsq(d$logSAM6, d$logCON6) #0.844
rsq(d$logSAM7, d$logCON7) #0.917
rsq(d$logSAM8, d$logCON8) #0.920
rsq(LFQ9$logSAM, LFQ9$logCON) #0.915
rsq(LFQ10$logSAM, LFQ10$logCON) #0.870

plot_grid(p1, p2, labels = c('C', 'D'), rel_widths = c(3, 2))
