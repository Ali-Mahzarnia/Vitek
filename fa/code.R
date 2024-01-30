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


 #noreadcsf=c(148,152,161,314,318,327) # dont read csf already in matlab

# 
# ### WHOLE BRAIN ANALYSIS
# # Specify and read folder with all of Nariman's data (individual_label_statistics/)
#  path_vol="/Users/ali/Desktop/Aug23/CVN/rba/exvivo/individual_label_statistics_ex_vivo/"
#  file_list=list.files(path_vol)

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

# Whole Brain Post-hoc Analysis
# Conduct post-hoc analysis and add to table AFTER ,0.05 and FDR correction
# Table 1
# Region -> Side of Brain -> Abbreviation -> index -> (3) Mean for each group -> (3) stedev for each group -> 
# FDR Corrected P-value -> F-value -> Eta-squared -> CI's
# Table 2
# Treatment -> T-ratio -> Uncorrected P-value -> Mean -> Standard Deviation -> Cohen Effect Size -> CI

######### create a similar format to voxelVolumes_treatments.xlsx for Volume:


# Create empty matrix to hold total brain Volume for each brain dataset
# temp=read.delim( paste0(path_vol,file_list[1]) )
len # =length(temp$Volume_mean)
 data_raw = matrix(NA,dim(BJsheet)[2], 334)

 
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
   data_raw[i,3:334]=tempppp / sum(tempppp)
   data_raw[i,2]= colnames(BJsheet)[(i+2)]
 }
 


# Set column names and filter vol_tab to only display CVN mice we are analyzing (by N-numbers)
struc_abbrev = read.csv('struc_abbrev.csv')



struc_temp = struc_abbrev
# struc_temp = struc_temp[-noreadcsf,]
struc_temp = struc_temp[1:332,]

colnames(data_raw) <- c('ID', 'N-number', struc_temp$Abbreviation)
data_raw <- data_raw[which(data_raw[,'N-number'] %in% metadata$Names),]
# data_raw2 <- data_raw2[which(data_raw2[,'N-number'] %in% metadata$N.number),]

# Replace N-numbers with each dataset's treatment condition for downstream ANOVA
for (g in 1:length(data_raw[,'N-number'])) {
  print(data_raw[,'N-number'][g])
  treat_row = which(data_raw[,'N-number'][g] == metadata$Names)
  data_raw[,'N-number'][g] = metadata$treatment[treat_row]
  print(metadata$treatment[treat_row])
  data_raw[,'ID'][g] = metadata$Names[treat_row]
}

# Rename 'N-number' column title with 'Treatment' column title; also convert Volume column
# to a numeric column type
colnames(data_raw)[2] = 'Treatment'
data_raw[,3:334] = as.numeric(data_raw[,3:334])
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


# remove the first 3 columns (index, filename, ID) and leave treatment and data columns; omit NaNs
data=na.omit(data) 
data_start = 5
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

# Remove first few columns
data=data[  , -c(1,2,3)]

## Create blank matrix to represent p values (column) for each brain region (rows)
colnames_vec = c("Full Structure", "Region ID", "Mean N ", "Mean Y", "SD N", "SD Y" 
                 ,"Uncorrected Pvalue", "FDR corrected Pvalue", "F-value", "Cohen's F Effect Size", "Effect Size Eta^2", 
                 "CI lower bound", "CI upper bound", 
                 "Shapiro-Wilk Pvalue (norm)", "Levene Test Pvalue (homog)")




pvalsresults=matrix(NA,(dim(data)[2]-1), length(colnames_vec))
rownames(pvalsresults)=names(data)[2:dim(data)[2]]
colnames(pvalsresults)=colnames_vec

# Populate created matrix with ANOVA-generated p-values
# Returning Partial Eta Squared (PES) effect size since we're using ANOVAs for each given brain region. Note
# that Cohen's d can be used however those are typically utilized when doing a t-test and evaluating a difference 
# between two groups
# Reference: https://www.youtube.com/watch?v=e5od9tH_QUo

# Comparison between ANOVA and linear model: 
# https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Learning_Statistics_with_R_-_A_tutorial_for_Psychology_Students_and_other_Beginners_(Navarro)/16%3A_Volumectorial_ANOVA/16.06%3A_ANOVA_As_a_Linear_Model
struc_abbrev = read.csv('struc_abbrev.csv')
len = dim(data)[2]
for (i in 1:(len-1))  {
  # Set 'tempname' to the currently analyzed brain regions
  tempname=rownames(pvalsresults)[i]
  
  # Find abbreviation full structure name
  sa_ind = which(struc_abbrev$Abbreviation == tempname)
  full_struc = struc_abbrev[sa_ind, 'Structure'][1]
  reg_ind = struc_abbrev[sa_ind, "Index"]

  # Create a linear model object using a given brain region and the data associated with the three
  # categories: 'sedentary', 'treadmill', 'wheel_only'. Calculate eta-squared effect size for the overall model 
  mylm <- lm(get(tempname) ~ as.factor(Treatment), data=data) 
  eff=effectsize::eta_squared(mylm, partial = F)
  cf=effectsize::cohens_f(mylm)
  
  # Conduct ANOVA on linear model to evaluate whether there is a significant difference between exercise treatment groups
  # within a given brain region; expect the means to be different and reject null (u_sedentary = u_treadmill = u_wheel-only)
  aov_table = anova(mylm)
  
  # Ho: data come from a normal distribution, H1: data do not come from a normal distribution
  # If p > 0.05, do NOT reject null, and thus data is normal
  normality = shapiro.test(mylm$residuals)
  
  # Ho: variances are equal, H1: at least one variance is different
  # If p > 0.05, do NOT reject null, and thus data has equal variances (homogeneity == equality of variances)
  homogeneity = leveneTest(get(tempname) ~ as.factor(Treatment), data=data)

  # Calculate the mean and standard deviation for each treatment group within the current brain region dataset
  data_temp = as.numeric( unlist(data[,i+1]))
  means=by(data_temp,as.factor(data$Treatment), mean)
  sds=by(data_temp,as.factor(data$Treatment), sd)
  
  # Output calculated values to the 'pvalresults' matrix in the user-defined order; the 'pvalresults' matrix will later
  # be used to create the 'posthoc_group_stats.csv' output
  val_list = c(full_struc, reg_ind, means[1], means[2], sds[1], sds[2], aov_table$`Pr(>F)`[1], aov_table$`Pr(>F)`[1], aov_table$'F value'[1], cf$Cohens_f, eff$Eta2, eff$CI_low, eff$CI_high, normality$p.value, homogeneity$`Pr(>F)`[1]) #normality$p.value>0.05, homogeneity$`Pr(>F)`[1]>0.05
  for (j in seq_along(val_list)){
    pvalsresults[i,j] <- val_list[j]
  }
}

colnames(pvalsresults)
pvalsresults[,8] = p.adjust(pvalsresults[,7], "fdr") #Benjamini & Hochberg

 write.csv( file = "results_volume.csv"  ,pvalsresults )
####


# Create new 'pvalsresultadjusted' variable to eventually populated with FDR-corrected values
pvalsresultsadjusted=pvalsresults

## adjust pvalues Benjamini & Hochberg
# To understand what we're doing here, we are comparing all the p-values determined from the ANOVAs conducted on all
# 332 brain regions. Because we can treat each ANOVA as an individual hypothesis test, we need to account for the
# multiple comparisons effect in the type 1 error rate. This can be done using various P-value (or significance)
# correction methods, but we will be using the Benjamini & Hochberg Correction method (specified by either 'BH or 'fdr')
# References: https://www.youtube.com/watch?v=rZKa4tW2NKs&t=483s
#             https://www.youtube.com/watch?v=K8LQSvtjcEo&t=43s
# pvalsresultsadjusted[,7] = p.adjust(pvalsresultsadjusted[,7], "fdr") #Benjamini & Hochberg
# pvalsresultsadjusted[,7] = as.numeric(pvalsresultsadjusted[,7])

# Filter 'pvalsresultesadjusted' table to display brain regions that have significant p values (p<0.05)
# sig = pvalsresults[as.numeric(pvalsresults[,10])<=0.05,]
sig = pvalsresults[pvalsresults[,7]<=0.05,] #&
                           # pvalsresultsadjusted[,ncol(pvalsresultsadjusted)]>0.05 & # For Levene's Test
                           # pvalsresultsadjusted[,ncol(pvalsresultsadjusted)-1]>0.05,] # For Shapiro-Wilk Test
sig = as.data.frame(sig)
# sig <- sig[order(sig$`Cohen's F Effect Size`, decreasing = TRUE),]


# Create empty 'posthoc' matrix to eventually populate with all gruop comparison data;
# Set first 'posthoc' column equal to 8th column of 'sig' matrix (FDR-corrected P-value column)
# posthoc = matrix(NA,dim(sig)[1],22)
# posthoc[,1]=sig[,"FDR corrected Pvalue"]
# 
# # Loop through each significant (FDR-corrected p-value < 0.05) brain region (from 'sig' matrix) and conduct
# # a Tukey test and report p-values for each comparison group
# for (i in 1:dim(sig)[1]) {
#   # Set 'tempname' to the currently analyzed brain regions
#   tempname=rownames(sig)[i]
# 
#   # Find abbreviation full structure name
#   sa_ind = which(struc_abbrev$Abbreviation == tempname)
#   full_struc = struc_abbrev[sa_ind, 'Structure']
#   reg_ind = struc_abbrev[sa_ind, "Index"]
#   #posthoc[,1 = full_struc]
# 
#   # Conduct ANOVA on currently analyzed brain region and use ANOVA output to conduct post-hoc Tukey test
#   mylm <- lm(get(tempname) ~ as.factor(Treatment), data=data)
#   # anova(mylm)
#   tuk=emmeans(mylm, list(pairwise ~ as.factor(Treatment)), adjust="tukey")
# 
#   # Get mean/stdev columns
# 
#   means = as.numeric(sig[tempname, c(3,4)])
#   sds = as.numeric(sig[tempname, c(5,6)])
# 
#   # Define Tukey comparison groups for calculating hedge's G effect size. Note that these groups are mannually
#   # set so that effect sizes can correlate with changes in regional FAs (e.g. negative effect size correlates
#   # with decrease in regional FAs)
#   control = c('N', 'Y')
#   treatment = c('N', 'Y')
# 
#   # For given brain region, loop through comparison groups and calculate hedge's g for each comparison groups
#   # for (j in 1:length(control)){
#   #   hedges_out = hedges_g(get(tempname) ~ factor(Treatment, levels=c(control[j], treatment[j])), data=data)
#   #   posthoc[i,j+13]=hedges_out$Hedges_g #Go through columns 5, 6, 7
#   # }
# 
#   # Add adjusted P-value and add low and high confidence intervals to 'posthoc' matrix
#   posthoc[i,2:4]=summary(tuk$`pairwise differences of Treatment`)$p.value
#   posthoc[i, 5:7]=summary(tuk$`pairwise differences of Treatment`)$t.ratio
#   posthoc[i, 8:10]=means
#   posthoc[i, 11:13]=sds
#   posthoc[i,c(17,19,21)] = summary(tuk$`emmeans of Treatment`)$lower.CL
#   posthoc[i,c(18,20,22)] = summary(tuk$`emmeans of Treatment`)$upper.CL
# }

# # Define output CSV column and row names
# colnames(posthoc)=c("FDR corrected Pvalue","SW Comparison Group Pvalue","ST Comparison Group Pvalue","TW Comparison Group Pvalue",
#                     "SW Comparison Group Tratio","ST Comparison Group Tratio","TW Comparison Group Tratio",
#                     "Mean sedentary group", "Mean voluntary group", "Mean voluntary + enforced group",
#                     "SD sedentary group", "SD voluntary group", "SD voluntary + enforced group",
#                     "SW Comparison Group Effect Size","ST Comparison Group Effect Size","TW Comparison Group Effect Size",
#                     "ST Comparison Group Lower CI", "ST Comparison Group Higher CI", "SW Comparison Group Lower CI", "SW Comparison Group Higher CI",
#                     "TW Comparison Group Lower CI", "TW Comparison Group Higher CI")
# rownames(posthoc)=rownames(sig)

# Save output CSV to specific folder string
# write.csv(sig, 'posthoc_group_stats_fa.csv')
# write.csv(posthoc, 'posthoc_comparison_stats_fa.csv')

# Activate Venv and run read_posthoc.py
# read_post.py serves to read above output CSV (posthoc_comparison_stats.csv specifically), and sort
# brain regions by effect sizes and postivitiy/negativity
data[,2:dim(data)[2]]=matrix_data_numeric/100
# 

# graph_list = list()

for ( jj in 1:dim(sig)[1]) {
  
          index= which(colnames(data) == rownames(sig)[jj] )
          p = ggplot(data, aes(x = Treatment, y = data[,index], fill = Treatment)) +
            geom_violin(alpha = 0.5) +
            geom_dotplot(binaxis = "y",
                         stackdir = "center",
                         dotsize = 0.5)  +
                 labs( y =  paste0(rownames(sig)[jj] ))
                ggsave( filename = paste0('plots/',rownames(sig)[jj],'.PNG' )  , p  )

}






# # Update Treatment names for plotting
# data_orig = data
# treat = unlist(data[,"Treatment"])
# # treat= gsub ("wheel_only", "Voluntary", unlist(treat ))
# # treat= gsub ("treadmill", "Voluntary + Enforced", treat )
# # treat= gsub ("sedentary", "Sedentary", treat )
# data[,"Treatment"] = treat
# # data[,"Treatment"] = data[,"Treatment"] %>% str_replace_all(c("wheel_only" = "Voluntary", "treadmill" = "Voluntary + Enforced", "sedentary" = "Sedentary"))
# 
# 
# # Define figure title hash/dict
# dict <- hash()
# dict[["NY"]] = c("N", "Y")
# # dict[["sw"]] = c("Sedentary", "Voluntary")
# # dict[["tw"]] = c("Voluntary", "Voluntary + Enforced")
# 
# ### PLOTTING
# # Set comparisong group vector to loop through: 'st': sedentary vs. treadmill, 'sw': sedentary vs. wheel_only
# # 'tw': treadmill vs. wheel_only
# comp_groups = c('NY')
# # plot_list = list()
# 
# 
# posthoc = as.data.frame(posthoc)
# library(formula.tools)

# Nested for loop: main loop specifies whether positive or negative effect sizes are being analysed
# inner for loop specifies which comparison group is being analyzed
# for (j in c( 'positive', 'negative')) {
  # print(paste(j, 'effect size regions'))
  # for (i in 1:length(comp_groups)) {
    # comparison = comp_groups[i] # 'st', 'sw', or 'tw'

    # Based on whether we are looking at the positive/negative effect size group, or the 'st', 'sw', or 'tw' comparison group,
    # we access the appropriate CSV in the below line; CSV is the output of the read_posthoc.py script

    #############
    # tempname1 = paste0( comparison, ' Comparison Group Pvalue')
    # ind1 = which(tempname1 == colnames(posthoc))
    # tempname2 = paste0( comparison, ' Comparison Group Effect Size')
    # ind2= which(tempname2 == colnames(posthoc))

    # if (j == 'positive') {       index = which(posthoc[,ind1]<0.05 & posthoc[,ind2]>0 )     }
    # if (j == 'negative') {       index = which(posthoc[,ind1]<0.05 & posthoc[,ind2]<0 )     }

    # if(length(index)>0){
    #   top_comp2 =  posthoc[index,]
      # # top_comp1 = matrix(NA, dim(top_comp2)[1], 10)
      # top_comp1[,1]= rownames(top_comp2)
      # top_comp[,1] =
      #




      #############
      # top_comp=read.csv(paste('/Users/ali/Desktop/Aug23/CVN/run_bass_our_fa/',comparison,'_regs/top_',substr(j,1,3),'_regions_',comparison,'.csv',sep=""))
      # data_comp = data[data$Treatment %in% dict[[i]],]
      # Assign values based on CSV columns to significant region names, abbreviations, and p-values and effect sizes
      # sig_reg = rownames(top_comp2)

      # indtemp = which(struc_temp$Abbreviation %in% sig_reg)
      # reg_struc = struc_temp$Structure[indtemp]
      # 
      # pvals_regs = formatC( as.numeric( top_comp2[,ind1], format = "e", digits = 2))
      # 
      # eff_sizes = formatC(as.numeric(top_comp2[,ind2], format = "e", digits = 2))

#       graph_list = list()
#       for (k in 1:3) {
#         if (is.na(sig_reg[k]) == F) {
#           res.aov <- aov(get(sig_reg[k]) ~ Treatment, data = data)
#           tuk=tukey_hsd(res.aov)
#           data_temp <- data[,c("Treatment",sig_reg[k])] %>% setNames(c("treatment","region"))
#           tuk <- add_y_position(tuk, data=data_temp, formula=region ~ treatment)
#           pbar_tab <- tuk[,c("group1", "group2", "p.adj", "y.position")]
# 
#           p <-  ggbetweenstats(
#             data= data,
#             x = Treatment,
#             y = !!sig_reg[k],
#             type = "parametric", # ANOVA or Kruskal-Wallis
#             var.equal = TRUE, # ANOVA or Welch ANOVA
#             plot.type = "box",
#             pairwise.comparisons = TRUE,
#             pairwise.display = "significant",
#             centrality.plotting = FALSE,
#             bf.message = FALSE,
#             ggsignif.args = list(textsize = 3, tip_length = 0.01),
#             smooth.line.args = list(linewidth = 4, color = "blue", method = "lm", formula = y ~x),
#             point.args = list(size = 4, alpha = 0.4, stroke = 0),
#             point.label.args = list(size = 4, max.overlaps = 1e+06),
#             ggplot.component = list(theme(text = element_text(size = 15))))
# 
# 
# 
# 
# 
#             # ggplot(data, aes_string(x="Treatment", y=sig_reg[k])) + stat_pvalue_manual(pbar_tab, label = "p.adj", size = 3, tip.length = 0, hide.ns = TRUE) +
#             # geom_violin() + geom_boxplot(width=0.1) + geom_dotplot(binaxis= "y", stackdir = "center", dotsize=0.5, fill='red') +
#             # labs(title=reg_struc[k], subtitle=paste("P-value of ",toString(pvals_regs[k])," | Effect size of ",toString(eff_sizes[k])), y="", x="") + theme_bw() +
#             # theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#             #
#           graph_list[[k]] = p
#         }
#         else {
#           p <- ggplot(data, aes_string(x="Treatment", y=sig_reg[length(sig_reg)])) + geom_blank() + theme_bw() + labs(x="", y="")
#           graph_list[[k]] = p
#         }
#       }
# 
#       p1 <- graph_list[[1]] #+ theme(axis.title.y = element_text(margin = margin(r = 20)), axis.title.x = element_blank())
#       p2 <- graph_list[[2]] #+ theme(axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_blank())
#       p3 <- graph_list[[3]] #+ theme(axis.title = element_blank())
# 
#       # Add total figure titles
#       patchwork <- p1 + p2 + p3 + plot_annotation(
#         #title = paste(dict[[comparison]][1],'vs.',dict[[comparison]][2],'Exercise'),
#         title = paste(comparison),
#         theme = theme(plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5, face = "bold"))
#       )
#       plot_list[[i]] = patchwork
#       full_plot <- p1 + p2 + p3 + plot_annotation(
#         title = paste('Regions of Significance following Post-Hoc Analysis (',str_to_title(j),' Effect Size)',sep=""),
#         # subtitle = paste(dict[[comparison]][1],'vs.',dict[[comparison]][2],'Exercise'),
#         subtitle = paste(comparison),
#         caption = 'DRAFT',
#         theme = theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#                       plot.margin=unit(c(1,1,-0.5,0), 'cm'))
#       )
#       gt <- patchwork::patchworkGrob(full_plot)
#       full_plot <- gridExtra::grid.arrange(gt, left = "  Volume (%)  ", bottom = "Treatment")
# 
#       # Saving individual plots
#       File <- paste(j,"_eff/",comparison,'_',substr(j,1,3),'.png',sep="")
#       ggsave(File, plot = full_plot, width=1213, height=514, dpi = 150, units='px', scale=2)
# 
#       # For tracking plotting and debugging
#       print(paste(comparison, 'complete'))
# 
#     }
#   }
#   if (length(plot_list) ==1){comp_plot = plot_grid(plot_list[[1]], nrow=1, ncol=1)}
#   if (length(plot_list) ==2){comp_plot = plot_grid(plot_list[[1]], plot_list[[2]], nrow=, ncol=1)}
#   if (length(plot_list) ==3){comp_plot = plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], nrow=3, ncol=1)}
#   composite_figure <- comp_plot + plot_annotation(
#     title = paste('Regions of Significance following Post-Hoc Analysis (',str_to_title(j),' Effect Size)',sep=""),
#     caption = 'DRAFT',
#     theme = theme(plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5),
#                   plot.margin=unit(c(1,1,-0.5,0), 'cm'))
#   )
#   gt <- patchwork::patchworkGrob(composite_figure)
#   composite_figure <- gridExtra::grid.arrange(gt, left = " Volume (%)", bottom = "Treatment")
#   # Saving Plots
#   File <- paste(substr(j,1,3),'.png',sep="")
#   ggsave(File, plot = composite_figure, width=1322, height=1322, dpi = 150, units='px', scale=2)
# }
# 
