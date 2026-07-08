## Author: Asia Jordan
### Date: 26/03/25
###: Female Urine: Regions - Total ! Not subsetting significant 
getwd()
setwd("/Users/asiajordan/Documents/2. PhD Bioinformatics/4. Urine Cohort/")

#Libraries 
library(ggfortify) 
library(minfi)
library(minfiData)
library(limma)
#install.packages("DMRcatedata_2.24.0.tar.gz", repos = NULL, type="source")
library(DMRcatedata)
library(DMRcate) #2.10.0 used for old epic, DMRcate_3.2.1.tar.gz for EPICv2
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(ggplot2)
library(beepr)
library(stargazer)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggdendro)
library(plyr)
library(dplyr)
library(RnBeads)
library(data.table)


getwd()
setwd("/Users/asiajordan/Dropbox/Asia PhD/4. Dissemination/Publications/Ch2. FemaleUrine/")

#bvals <- read.csv("CarsonAnalysis/menopause_urine_out/betas_full.csv")
#colnames(bvals) <- c("IlmnID","FMLS_007", "FMLS_009", "FMLS_005","FMLS_006","FMLS_011","FMLS_015", "FMLS_002","FMLS_010")
#rownames(bvals) <- bvals$IlmnID
#bvals$IlmnID <- NULL
#head(bvals)

#originalmval <- read.csv("/Users/asiajordan/Documents/2. PhD Bioinformatics/4. Urine Cohort/Female_Urine_Methylation_Analysis_Dermot_multiComp/Results/mVals_AllCpGs_Urine.csv")
mval <- read.csv("CarsonAnalysis/scripts/Asia/6_CarsonMVals_ALL.csv", row.names = 1)

 #beep()

#Phenotypic data
ClinicData = read.csv("CarsonAnalysis/EPIC_FemaleUrine_SampleSheet.csv")

#Load comparison data
Urine <- read.csv("CarsonAnalysis/scripts/Asia/2_DiffMeth_DMPs_EPICdata_BetaValDif_Significant.csv")
nrow(Urine) #149649
colnames(Urine)
#beep()
filtered_data_dif <- Urine %>%
  filter((BetaValueDif >= 0.2 | BetaValueDif <= -0.2) & adj.P.Val <=0.05)
dim(filtered_data_dif) # 149649     21

#Get significant probes
DMPs = filter(Urine, adj.P.Val<0.05 & BetaValueDif>=0.2 | BetaValueDif<=-0.2)
dim(DMPs) #149649     21
head(DMPs)
DMPs_Sig <- as.data.frame(DMPs[,2])
rownames(DMPs_Sig) <- DMPs_Sig$`DMPs[, 2]`

# get mvalues for significant DMPs
DMPs_Sig_mvals <- merge(DMPs_Sig, mval, by=0)
rownames(DMPs_Sig_mvals) <- DMPs_Sig_mvals$Row.names
dim(DMPs_Sig_mvals) # 149649     10
head(DMPs_Sig_mvals)
DMPs_Sig_mvals <- DMPs_Sig_mvals[,3:10] #Need all numeric values
# use the ClinicData dataframe to create a design matrix

#Match Order
ClinicData <- ClinicData[
  match(colnames(mval), ClinicData$Sample.ID),
]
all(ClinicData$Sample.ID == colnames(mval))
ClinicData$Sample.ID
colnames(mval)
ClinicData
head(mval)

ClinicData <- ClinicData %>%
  mutate(Menopause.stage = recode(Menopause.stage, 
                                  "Post-menopause" = "PostMeno", 
                                  "Pre menopause" = "PreMeno"))
Targets <- unique(ClinicData$Menopause.stage) #Comparison options
Targets
Group <- factor(ClinicData$Menopause.stage, levels=Targets) #Groups
Group
Design <- model.matrix(~0+Group) #Design matrix
Design
colnames(Design) <- c("PostMeno", "PreMeno") #Rename columns of design matrix

#Create contMatrix (statistical linear model fit using M values and design matrix (adjusted p values are calculated))
fit <- lmFit(mval, Design) #Model
fit
#beep()
contMatrix <- makeContrasts(PostMeno-PreMeno, levels=Design)  #contMatrix 
contMatrix
#Create annEPIC manifesto for annotation

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
?getAnnotation()
class(annEPIC)
dim(annEPIC) #930075     42
colnames(annEPIC)
annEPICSub <- annEPIC[match(rownames(mval),annEPIC$Name), c(1:4,9,18:19,21,25:28,31, 34, 38:ncol(annEPIC))] #Pulling out the columns we want in final dataframes - NB: UCSC gene details, Fantom4 and Fantom 5, 450K array, 
colnames(annEPICSub)
head(annEPICSub)

#### DMRcate ####
#Looking at all probes
?cpg.annotate()
mval[1:9,1:8]
?rmSNPandCH()
dim(mval) # 886928      8
type(mval)
class(mval)
mval <- as.matrix(mval)
head(mval)

myAnnotation1 <- cpg.annotate(object = as.matrix(DMPs_Sig_mvals), 
                              datatype = "array", 
                              what = "M",
                              analysis.type = "differential", 
                              design = Design,
                              epicv2Remap = FALSE,
                              epicv2Filter = "precision", #Options: Strategy for filtering probe replicates that map to the same CpG site. "mean" takes the mean of the available probes; "sensitivity" takes the available probe most sensitive to methylation change; "precision" either selects the available probe with the lowest variation from the consensus value (most precise), or takes the mean if that confers the lowest variation instead, "random" takes a single probe at random from each replicate group.
                              contrasts = TRUE, 
                              cont.matrix = contMatrix,
                              coef = "PostMeno - PreMeno",  # nad to put a " " either side of the '-' to fix error: Error in cpg.annotate(object = as.matrix(mval), datatype = "array", what = "M",  : coef %in% colnames(cont.matrix) is not TRUE
                              arraytype = "EPICv2", 
                              #epicv2Filter collapses replicate probes by taking their mean.
                              #object=ALLMs.noSNPs,
                              fdr = 0.05)
#Your contrast returned 148216 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
#beep()
#Your contrast returned 115487 individually significant probes. We recommend the default setting of pcutoff in dmrcate().

DMRs_sigSite <- dmrcate(myAnnotation1, lambda=1000, C=2, min.cpgs=3 )  #min.cpgs=Minimum number of consecutive CpGs constituting a DMR
str(DMRs_sigSite)
DMRs_sigSite_ranges = extractRanges(DMRs_sigSite) 
#beep()
DMRs_sigSite_list = as.data.frame(DMRs_sigSite_ranges)
dim(DMRs_sigSite_list) #8881   13
head(DMRs_sigSite_list)

#### Extract CpG sites from DMRs and annotate results from DMRcate ####
#Total number of CpGs that were selected for DMRs
#Step 1: Make the dataframe
dim(annEPICSub)
rownames(annEPICSub) <- make.unique(
  ifelse(is.na(rownames(annEPICSub)),
         "missing",
         rownames(annEPICSub))
)
annEPICSub2 <- as.data.frame(annEPICSub)
DMRs_sites_list <- data.frame(matrix(nrow = 0, ncol = ncol(annEPICSub2)))
colnames(DMRs_sites_list)=colnames(annEPICSub2)
dim(annEPICSub2)

head(DMRs_sigSite_list) #seqnames, start, end, 
head(annEPICSub2) # chr, pos, Name, 
head(DMRs_sites_list) #Matches above

#Loop to find DMPs for each DMR
for(i in 1:nrow(DMRs_sigSite_list)){ #For loop iterating over all rows in DMRs_sigSite_list
  Resultscpg <- annEPICSub2 %>%  #The next part subsets rows from annEPICSub2 based on the following:
    subset(chr == DMRs_sigSite_list$seqnames[i] & pos %in% DMRs_sigSite_list$start[i]:DMRs_sigSite_list$end[i]) %>% #1) value of chr column matches annEPICSub2 seqnames column, 2) pos (position) values of annEPICSub2 are within the range defined by DMRs_sigSite_list$start and DMRs_sigSite_list$end  
    arrange(pos) #After the above, this sorts the dataframe by pos (position) in ascending order
  DMRs_sites_list <- rbind(DMRs_sites_list, Resultscpg) #Then in each iteration of the loop the Resultscpg dataframe is appended to the DMRs_sites_list dataframe 
}
#beep()
head(DMRs_sigSite_list)
head(DMRs_sites_list) #Take a look at what the loop made


#Extract cpgs from the end region and collapse rows into one for annotation
Results <- data.frame(matrix(nrow = nrow(DMRs_sigSite_list), ncol = ncol(annEPICSub2)))
colnames(Results) = colnames(annEPICSub2)

for(i in 1:nrow(DMRs_sigSite_list)){
  Results[i,] <- annEPICSub2 %>% 
    subset(chr == DMRs_sigSite_list$seqnames[i] & pos %in% DMRs_sigSite_list$start[i]:DMRs_sigSite_list$end[i]) %>%
    arrange(pos) %>% #The start of this loop is the same as the above but instead creates results dataframe (not resesultscpg)
    lapply(paste0, collapse=",") #And this time instead of above, we are applying the paste0 function to concatenate elements together without spaces, and it's collapsing them into a single string separated by commas. 
}  
#beep()
head(DMRs_sigSite_list) 
head(Results)  #Take a look at what the loop made

Results <- Results[,-c(1:2)] #Remove repetitive chromosome no and positional info
DMRs_sigSite_list <- tibble::rowid_to_column(DMRs_sigSite_list, "n") #Add new column with unique identifier names n
Results <- tibble::rowid_to_column(Results, "n")
head(Results)
DMRs_subset <- join(DMRs_sigSite_list, Results, by="n") #Merge by this new column we created
DMRs_subset <- DMRs_subset[,-1] #Remove n column

#Assign a DMR ID 
DMRs_subset=data.frame(paste("MENO_DMR",c(1:nrow(DMRs_subset)),sep="_"), DMRs_subset)
colnames(DMRs_subset)[1] = "DMR_ID"

#select only significant DMRs & write these to a csv
colnames(DMRs_subset) 
DMRs_subset = filter(DMRs_subset, HMFDR<0.05 & meandiff>=0.2 | meandiff<=-0.2 )
dim(DMRs_subset)  #A subset of DMRs are significant and different >0.1
# 8213   31
write.csv(DMRs_subset,"CarsonAnalysis/scripts/Asia/7_1_DMRcate_EPICv2_Menopause_Urine_SignificantDMRs.csv") 

#Uncollapse DMRs into DMP IDs
DMRs = DMRs_subset[,c("DMR_ID", "Name")] #Pull out two columns from DMR_subset
DMP_DMR <-  DMRs %>%
  separate_longer_delim(c(1:ncol(DMRs)), delim = ",")
head(DMP_DMR)
table(duplicated(DMP_DMR$DMR_ID)) # unique 6102 duplicated 23046, Note 23046 needs to be divided to get the actual unique number of CpGs
DMP_DMR <- as.data.frame(DMP_DMR)
DMPs <- as.data.frame(DMPs)
dim(DMP_DMR) #38943     2
#Write csv for list of CpG probe IDs
write.csv(DMP_DMR,"CarsonAnalysis/scripts/Asia/7_2_DMRcate_EPICv2_Menopause_Urine_SignificantDMRs_DMPs.csv")  

#For a more details on each CpG 
rownames(DMP_DMR) <- DMP_DMR$Name
rownames(DMPs) <- DMPs$IlmnID
DMPs = merge(DMPs, DMP_DMR, by=0)
dim(DMPs) # 34439    24
write.csv(DMPs,"CarsonAnalysis/scripts/Asia/7_3_DMRcate_EPICv2_Menopause_Urine_SignificantDMRs_DMPs_Details.csv") 

#For a far more in depth overview of each CpG, this is from details from the EPIC manifest B4 annotation (hg19)
DMPs <- read.csv("CarsonAnalysis/scripts/Asia/7_3_DMRcate_EPICv2_Menopause_Urine_SignificantDMRs_DMPs_Details.csv")
rownames(DMPs) <- DMPs$IlmnID
DMPs$cgid <- NULL
head(annEPICSub2)
DMPsMan = merge(DMPs, annEPICSub2, by=0)
dim(DMPsMan) #34439    45
EPIC_edited <- read.csv("/Users/asiajordan/Documents/2. PhD Bioinformatics/4. Urine Cohort/MethylationEPIC v2.0 Files/EPICv2_Annotation_Summerised_AllProbeTypes_EnhFixed.csv")
EPIC_edited[0:3,0:3]
table(EPIC_edited$Enhancer)
rownames(EPIC_edited) <- EPIC_edited$IlmnID
EPIC_edited$IlmnID<- NULL
DMPsMan[0:3,0:3]
rownames(DMPsMan) <- DMPsMan$Row.names
DMPsMan$Row.names <- NULL
DMPsMan1 = merge(DMPsMan, EPIC_edited, by=0)
dim(DMPsMan1) # 34439    82
table(DMPsMan1$Enhancer)
write.csv(DMPsMan1,"CarsonAnalysis/scripts/Asia/7_4_DMRcate_EPICv2_Menopause_Urine_SignificantDMRs_DMPs_MANIFESTDetails_EnhFixed.csv") 

