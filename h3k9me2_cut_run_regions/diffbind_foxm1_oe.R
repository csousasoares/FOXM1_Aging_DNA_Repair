## Load packages and set seed

set.seed(123)

library(DiffBind)
library(EnrichedHeatmap)
library(colorRamp2)
library(tidyverse)
library(profileplyr)

## Load sample info

samples <- read.csv("h3k9me2.csv", sep = ";")
names(samples)
samples

## Create a DBA object

h3k9me2 <- dba(sampleSheet = samples)
h3k9me2

plot(h3k9me2)

## Count reads in re-centered peaks/regions

h3k9me2 <- dba.count(h3k9me2, summits = 1000, minOverlap = 2) 

## Since H3K9me2 spans broad regions, regions were analysed in
## large 2Kb (2 x 1000 bp) intervals
## Set minOverlap = 2 to make sure regions are reproducible in at least 
## two samples (most likely replicates of either condition)

h3k9me2
plot(h3k9me2)

#Normalize the data

## First normalize by sequencing depth, as if there were no differences
## between samples. This is important because sequencing depth alters the 
## DiffBind estimate of Log2FC, even if samples are equivalent.

## When all samples have the same signal, obtain the true signal by using
## the spike-in calibration ratio

spikeinreads <- c(1,1,1,1) 

## First, do 1,1,1,1, to normalize by depth.
## Then, replace by normalization_method_final

length(spikeinreads)

h3k9me2 <- dba.normalize(h3k9me2, normalize = spikeinreads)
norm_spikes <- dba.normalize(h3k9me2, normalize = spikeinreads, bRetrieve = T)
norm_spikes

normlibs_spike <- cbind(FullLibSize = norm_spikes$lib.sizes, 
                        NormFacs = norm_spikes$norm.factors,
                        NormLibSize = round(norm_spikes$lib.sizes/norm_spikes$norm.factors))
normlibs_spike

df_normlibs_spike <- as.data.frame(normlibs_spike)

normalizatin_method <- df_normlibs_spike$NormLibSize/max(df_normlibs_spike$NormLibSize)
normalizatin_method

## Multiply this by spike-in ratio:

spike_in_ratios <- c(1,1,0.42,0.63)

normalizatin_method_final <- normalizatin_method*spike_in_ratios
normalizatin_method_final

spikeinreads <- c(1, 0.493, 0.415, 0.37)
length(spikeinreads)

h3k9me2 <- dba.normalize(h3k9me2, normalize=spikeinreads)
norm_spikes <- dba.normalize(h3k9me2, normalize=spikeinreads, bRetrieve = T)
norm_spikes

normlibs_spike <- cbind(FullLibSize=norm_spikes$lib.sizes, 
                        NormFacs=norm_spikes$norm.factors,
                        NormLibSize = 
                          round(norm_spikes$lib.sizes/norm_spikes$norm.factors))
normlibs_spike

df_normlibs_spike <- as.data.frame(normlibs_spike)

normalizatin_method <- df_normlibs_spike$NormLibSize/
  max(df_normlibs_spike$NormLibSize)


## Now treated samples have aprox. the same sequencing depth
## while the sequencing depth of the control samples
## is spike-in corrected relative to the treated depth

## This is done relative to treated samples to avoid artificially
## reducing the number of reads to decimal values between 0 and 1.
## It also prevents Diffbind from "thinking" that control samples have
## different binding intensity from one another.


## Set the contrasts like in DESeq2:

h3k9me2$contrasts <- NULL ## Ensures it is empty

h3k9me2 <- dba.contrast(h3k9me2, categories = DBA_TREATMENT, minMembers=2, 
                        contrast = c("Treatment", "DNDK", "VECTOR"))
h3k9me2

##Analyze differential binding between conditions with DESeq2

h3k9me2_analysis <- dba.analyze(h3k9me2, method = DBA_DESEQ2, 
                                bGreylist = F) 

## DESeq2 method works nicely
## Blacklist is applied by default


dba.plotMA(h3k9me2_analysis) ## Show MA Plot

dba.show(h3k9me2_analysis, bContrasts=TRUE) ## Check the contrasts

h3k9me2.DB <- dba.report(h3k9me2_analysis, th = 1) 
##Get GRanges - will be transformed in df later
h3k9me2.DB

dba.plotBox(h3k9me2_analysis) ## BoxPlot of Log2(reads) in peaks
dba.plotVolcano(h3k9me2_analysis) ## Volcano Plot of Regions

##Create dataframe with results

final_results <- as.data.frame(h3k9me2.DB)

final_results$bh <- p.adjust(final_results$p.value, method = "BH")

write.csv(final_results, "awesome_results_foxm1_oe.csv")

##Creates PCA

dba.plotPCA(h3k9me2_analysis,DBA_TREATMENT,label=DBA_TREATMENT)

# Capture the plot object
p <- dba.plotPCA(h3k9me2_analysis, DBA_TREATMENT, label=DBA_TREATMENT)

# Modify the aspect ratio
p$aspect.ratio <- 0.5

# Plot it
print(p)

hmap <- colorRampPalette(c("deepskyblue2", "white", "red"))(n = 13)


readscores <- dba.plotHeatmap(h3k9me2_analysis, contrast=1, 
                              correlations=FALSE, scale="row", 
                              colScheme = hmap, th = 0.05)



#----------------------#
#-----Plot Profile-----#
#----------------------#

?dba.plotProfile

repObj <- dba.report(h3k9me2_analysis, contrast = 1, bDB=TRUE)
profiles <- dba.plotProfile(h3k9me2_analysis, sites = repObj, 
                            maxSites = 1000, 
                            window_size = 2000, 
                            style = "percentOfRegion", 
                            distanceAround=30)


cols <- colorRamp2(c(0, 3, 6, 10), viridis::viridis(4))
args(dba.plotProfile)
profiles
rowRanges(profiles)[1:2]
heatmap <- generateEnrichedHeatmap(profiles, raster_quality = 15, 
                                   ylim = c(0, 12), 
                                   matrices_color = cols, 
                                   use_raster = length(profiles) > 500)
dba.plotProfile
