### LIBRARY IMPORTS
library(data.table)
library(readxl)
library(rstatix)
library(jmv)
library(ggplot2)
library(patchwork)
library(ggeasy)
library(hash)
library(emmeans)
library(effectsize)
library(reticulate)
library(carData)
library(car)
library(stringr)
library(cowplot)
library(ggpubr)
library(dplyr)
library(ggstatsplot)


# Specify and read ID_Treatment.CSV - this CSV associates contains metadata for all the mice brain
# data files: Original ID, Modified ID, N-number, Treatment
path_metadata = "VITEK_metadata.xlsx"
metadata = read_xlsx(path_metadata, sheet = 2)

# *******
# Remove 'A19120503_T1' ID dataset since we don't have that file
# metadata = metadata[!(metadata$Modified.ID=='A19120503_T1'),]
# *******

# Create empty matrix to hold total brain FA for each brain dataset
# temp=read.delim( paste0(path_vol,file_list[1]) )


## BJ sheet
BJsheet  = read.delim('studywide_stats_for_volume.txt')
dim(BJsheet)

len=dim(BJsheet)[1]
vol_tab = matrix(NA,(dim(BJsheet)[2]-2),3)

# Populate matrix with brain dataset filename and total FA (mm3)
for (i in 1:(dim(BJsheet)[2]-2)) {
  # temp=read.delim(paste0(path_vol,file_list[i]))
  vol_tab[i,3]=sum( as.numeric(BJsheet[2:len,(i+2)])  )
  vol_tab[i,2] = colnames(BJsheet)[(i+2)]
  
  
}


# Set column names and filter vol_tab to only display CVN mice we are analyzing (by N-numbers)
colnames(vol_tab) <- c('ID', 'N-number', 'Volume')
vol_tab <- vol_tab[which(vol_tab[,'N-number'] %in% metadata$Names),]

# Replace N-numbers with each dataset's treatment condition for downstream ANOVA
for (g in 1:length(vol_tab[,'N-number'])) {
  treat_row = which(vol_tab[,'N-number'][g] == metadata$Names)
  vol_tab[,'N-number'][g] = metadata$treatment[treat_row]
  vol_tab[,'ID'][g] = metadata$Names[treat_row]
}

# Rename 'N-number' column title with 'Treatment' column title; also convert Volume column
# to a numeric column type
colnames(vol_tab)[2] = 'Treatment'
wb_data = data.frame(vol_tab)
wb_data[,3] = as.numeric(wb_data[,3])

# Conduct and report whole-brain ANOVA
wblm = lm(Volume ~ Treatment, data=wb_data[,c('Treatment','Volume')])
wb_aov = anova(wblm)

# Print TRUE if anova p-value < 0.05
paste('Whole-brain_Volume:',(wb_aov$`Pr(>F)` < 0.05)[1])

# Save Whole-brain ANOVA data 
wb_data = wb_data[order(wb_data$Treatment, decreasing = TRUE), ]
write.csv(wb_data, 'wb_Volume.csv', row.names = F)



# Create empty matrix to hold total brain Volume for each brain dataset
# temp=read.delim( paste0(path_vol,file_list[1]) )
len # =length(temp$Volume_mean)
data_raw = matrix(NA,dim(BJsheet)[2], 336)


# Populate matrix with brain dataset filename and total Volume (mm3)
#  for (i in 1:(dim(BJsheet)[2]-2)) {
#   # temp=read.delim(paste0(path_vol,file_list[i]))
#   tempppp =  as.numeric(BJsheet[2:len,(i+2)])  
#   data_raw[i,3:334]=tempppp
#   data_raw[i,2]=substr(file_list[i], 1, 6)
# }

# data_raw = matrix(NA,dim(BJsheet)[2], 334)

# Populate matrix with brain dataset filename and total Volume (mm3)
for (i in 1:(dim(BJsheet)[2]-2)) {
  # temp=read.delim(paste0(path_vol,file_list[i]))
  tempppp =  as.numeric(BJsheet[2:len,(i+2)])  
  data_raw[i,5:336]=tempppp / sum(tempppp)
  data_raw[i,2]= colnames(BJsheet)[(i+2)]
}



# Set column names and filter vol_tab to only display CVN mice we are analyzing (by N-numbers)
struc_abbrev = read.csv('struc_abbrev.csv')



struc_temp = struc_abbrev
# struc_temp = struc_temp[-noreadcsf,]
struc_temp = struc_temp[1:332,]

colnames(data_raw) <- c('ID', 'N-number' ,  'Genotype', 'Sex', struc_temp$Abbreviation)
data_raw <- data_raw[which(data_raw[,'N-number'] %in% metadata$Names),]
# data_raw2 <- data_raw2[which(data_raw2[,'N-number'] %in% metadata$N.number),]

# Replace N-numbers with each dataset's treatment condition for downstream ANOVA
for (g in 1:length(data_raw[,'N-number'])) {
  print(data_raw[,'N-number'][g])
  treat_row = which(data_raw[,'N-number'][g] == metadata$Names)
  data_raw[,'N-number'][g] = metadata$treatment[treat_row]
  data_raw[,'Genotype'][g] = metadata$Genotype[treat_row]
  data_raw[,'Sex'][g] = metadata$Sex[treat_row]
  
  print(metadata$treatment[treat_row])
  data_raw[,'ID'][g] = metadata$Names[treat_row]
}

# Rename 'N-number' column title with 'Treatment' column title; also convert Volume column
# to a numeric column type
colnames(data_raw)[2] = 'Treatment'
data_raw[,5:336] = as.numeric(data_raw[,5:336])
data_raw = data.frame(data_raw)
data_raw = as.data.frame(cbind(data_raw$ID ,   data_raw$ID, data_raw  ))
colnames(data_raw)[1] = 'Entry'
colnames(data_raw)[2] = 'Filename'
data=data_raw

### BEGIN STATISTICAL ANALYSIS
# Read voxel Volumes w/ treatments excel file
# data=read_xlsx('/Users/ali/Desktop/Aug23/CVN/Bass-Connections-F22-main/Statistical_Analysis/voxelVolumes_treatments.xlsx')

# Replace appended "X" to region names due to read.csv wrapper
old_colnames = colnames(data)[substr(colnames(data),1,1)=="X"]
new_colnames = sub('.','',old_colnames)
colnames(data)[colnames(data) %in% old_colnames] <- new_colnames

# 

# remove the first 3 columns (index, filename, ID) and leave treatment and data columns; omit NaNs
data=na.omit(data) 
data_start = 7
dim(data)

matrix_data_numeric = data_raw[,data_start:dim(data)[2]]
matrix_data_numeric = matrix_data_numeric %>%  mutate_all(as.numeric)
# Convert voxel counts to proportions of total brain Volume
data[,data_start:dim(data)[2]]=100*matrix_data_numeric




# Add Whole-brain data to conduct full analysis
wb_vec = matrix(NA,length(data$ID),2)
for (i in 1:length(data$ID)) {
  val = which(data$ID[i] == vol_tab)
  wb_vec[i, 1] = data$ID[i]
  wb_vec[i, 2] = vol_tab[val, 'Volume']
}
colnames(wb_vec) = c('ID', 'Volume')
data$Brain <- as.numeric(wb_vec[,'Volume'])
dim(data)


pvalsresults=matrix(NA,(dim(data)[2]-3), 20)
rownames(pvalsresults)=names(data)[4:dim(data)[2]]
# colnames(pvalsresults)=colnames_vec

# Populate created matrix with ANOVA-generated p-values
# Returning Partial Eta Squared (PES) effect size since we're using ANOVAs for each given brain region. Note
# that Cohen's d can be used however those are typically utilized when doing a t-test and evaluating a difference 
# between two groups
# Reference: https://www.youtube.com/watch?v=e5od9tH_QUo

# Comparison between ANOVA and linear model: 
# https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Learning_Statistics_with_R_-_A_tutorial_for_Psychology_Students_and_other_Beginners_(Navarro)/16%3A_Volumectorial_ANOVA/16.06%3A_ANOVA_As_a_Linear_Model
struc_abbrev = read.csv('struc_abbrev.csv')
pvalsresults = matrix(NA,(dim(data)[2]), 20)
for (i in 7:dim(data)[2]) {
  
  tempname=colnames(data)[i]
  
  # Find abbreviation full structure name
  sa_ind = which(struc_abbrev$Abbreviation == tempname)
  full_struc = struc_abbrev[sa_ind, 'Structure'][1]
  reg_ind = struc_abbrev[sa_ind, "Index"][1]
  
  # Create a linear model object using a given brain region and the data associated with the three
  # categories: 'sedentary', 'treadmill', 'wheel_only'. Calculate eta-squared effect size for the overall model 
  mylm <- lm(get(tempname) ~ as.factor(Genotype)*as.factor(Treatment)*as.factor(Sex), data=data) 
  # eff=effectsize::eta_squared(mylm, partial = F)
  # cf=effectsize::cohens_f(mylm)
  aov_table = anova(mylm)
  
  normality = shapiro.test(mylm$residuals)
  homogeneity = leveneTest(get(tempname) ~ as.factor(Treatment), data=data)
  data_temp = as.numeric( unlist(data[,i]))
  val_list = c(full_struc, reg_ind, aov_table$`Pr(>F)`, aov_table$'F value', normality$p.value, homogeneity$`Pr(>F)`[1]) #normality$p.value>0.05, homogeneity$`Pr(>F)`[1]>0.05
  colnames(pvalsresults) = c("Region","Region ID" ,paste0("P-value ", rownames(aov_table)),rownames(aov_table), "P-value Normality", "P-value Homogeneity"   )
  for (j in seq_along(val_list)){
    pvalsresults[i,j] <- val_list[j]
  }
  
}
pvalsresults = pvalsresults[-seq(1,6),] 
pvalsresultsadjusted = pvalsresults

index_pval_correction = grep( "P-value", colnames( pvalsresults))

for (j in index_pval_correction) {
  pvalsresultsadjusted[,j] = p.adjust(pvalsresultsadjusted[,j], "fdr") #Benjamini & Hochberg
}


write.csv( file = "results_volume_ANOVA.csv"  ,pvalsresultsadjusted )




