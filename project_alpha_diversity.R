##############################################################
# title: "Alpha diversity in R - qiime2 output"
# author: "ANSC595"
# date: "March 16, 2021"
##############################################################

#The first step is very important. You need to set your working 
#directory. Just like in unix we have to `cd` to where our data is, the 
#same thing is true for R.
##############################################

setwd("F:/ANSC516/PROJECT/")
list.files()

# Modified from the original online version available at 
# http://rpubs.com/dillmcfarlan/R_microbiotaSOP

# and Tutorial: Integrating QIIME2 and R for data visualization 
# and analysis using qiime2R
# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121



##Files
# We will use the following files created using qiime2

# evenness_vector.qza (alpha diversity)
# faith_pd_vector.qza (alpha diversity)
# observed_features_vector.qza (alpha diversity)
# shannon_vector.qza (alpha diversity)
# MetaData_complete2.txt (metadata)


# Data manipulation
## Load Packages

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","

# Load data

meta<-read_q2metadata("MetaData_complete2.txt")
str(meta)
#colnames(meta)[3] <- "new name for the column"
#str(meta)

evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd = faith_pd[,-1]
colnames(faith_pd)[1] <- "SampleID"
colnames(faith_pd)[2] <- "faith_pd"


## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(meta)

str(observed_features)



###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta_merged = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta_merged) <- meta_merged$SampleID #this makes a first column on row names and 
#names each of the rows the sample id associated with it 
meta_merged = meta_merged[,-1] #removes the sample.id column so your table is less redundant (this isn't necessary at all)
str(meta_merged)


#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or use non-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(meta_merged$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta_merged$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta_merged$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta_merged$observed_features), main="Observed Features", xlab="", breaks=10)
############none appear normal

#Plots the qq-plot for residuals
ggqqplot(meta_merged$shannon_entropy, title = "Shannon")
ggqqplot(meta_merged$faith_pd, title = "Faith PD")
ggqqplot(meta_merged$pielou_e, title = "Evenness")
ggqqplot(meta_merged$observed_features, title = "Observed Features")

#install.packages("ggpubr")
library("ggpubr")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta_merged$shannon)
shapiro.test(meta_merged$faith_pd)
shapiro.test(meta_merged$pielou_e)
shapiro.test(meta_merged$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest:

## Categorical variables
# Now that we know which tests can be used, let's run them. 

## Normally distributed metrics

# for a categorical variable with more than 2 levels, we run ANOVA. 
#If it's only two levels (binary), we could run a t-test

#####################not using this for project, as none are normal################
#Run the ANOVA and save it as an object
#aov.evenness.body_site = aov(pielou_evenness ~ body.site, data=meta)
#Call for the summary of that ANOVA, which will include P-values
#summary(aov.evenness.body_site)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

#TukeyHSD(aov.evenness.body_site)

# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.

#levels(meta$body.site)
#Re-order the groups because the default is alphabetical order
#you have to reorder it now so that it comes out in the right order on the graph
#eta$body.site.ord = factor(meta$body.site, c("left palm", "right palm", "gut", "tongue"))
#levels(meta$body.site.ord)

#Plot
#boxplot(pielou_evenness ~ body.site.ord, data=meta, ylab="Peilou evenness")

#evenness_boxplot <- ggplot(meta, aes(body.site.ord, pielou_evenness)) + # _boxplot uses median and percentile, _se uses avg and standard error
  #geom_boxplot() + 
  #ylim(c(0.5,1)) +
  #theme_q2r() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#ggsave("output/evenness_boxplot.png", evenness_boxplot, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. 
#ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. 
#Unfortunately plotting the average and standard deviation is a little complicated.

#evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  #group_by(body.site.ord) %>%   # the grouping variable
  #summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            #sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            #n_evenness = n(),  # calculates the sample size per group
            #se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

#evenness_se <- ggplot(evenness_summary, aes(body.site.ord, mean_evenness, fill = body.site.ord)) + 
  #geom_col() + 
  #geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  #theme_q2r() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #theme(legend.title = element_blank()) +
  #labs(y="Pielou's evenness  ± s.e.", x = "") 

#ggsave("output/evenness_se.png", evenness_se, height = 2.5, width = 3)
################################################################################################

## **Non-normally distributed metrics**

# We will use **Faith's phylogenetic diversity** here. 
# for categorical variables, we use Kruskal-Wallis (non-parametric equivalent of ANOVA). 
# If we have only two levels, we would run Wilcoxon rank sum test (non-parametric equivalent of t-test)

kruskal.test(faith_pd ~ Tissue, data=meta_merged)
kruskal.test(faith_pd ~ AgeAtHarvest, data=meta_merged)
kruskal.test(pielou_evenness ~ Tissue, data=meta_merged)
kruskal.test(pielou_evenness ~ AgeAtHarvest, data=meta_merged)
kruskal.test(observed_features ~ Tissue, data=meta_merged)
kruskal.test(observed_features ~ AgeAtHarvest, data=meta_merged)
kruskal.test(shannon_entropy ~ Tissue, data=meta_merged)
kruskal.test(shannon_entropy ~ AgeAtHarvest, data=meta_merged)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(meta_merged$faith_pd, meta_merged$`EcoliTreatment`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$faith_pd, meta_merged$`SeedSanitization`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$pielou_evenness, meta_merged$`EcoliTreatment`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$pielou_evenness, meta_merged$`SeedSanitization`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$observed_features, meta_merged$`EcoliTreatment`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$observed_features, meta_merged$`SeedSanitization`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$shannon_entropy, meta_merged$`EcoliTreatment`, p.adjust.method="BH")
pairwise.wilcox.test(meta_merged$shannon_entropy, meta_merged$`SeedSanitization`, p.adjust.method="BH")

#Plot
boxplot(faith_pd ~ Tissue, data=meta_merged, ylab="Faith phylogenetic diversity")
boxplot(faith_pd ~ EcoliTreatment, data=meta_merged, ylab="Faith phylogenetic diversity") ###########not working????
boxplot(faith_pd ~ SeedSanitization, data=meta_merged, ylab="Faith phylogenetic diversity")
boxplot(faith_pd ~ AgeAtHarvest, data=meta_merged, ylab="Faith phylogenetic diversity")

# or with ggplot2
tissue_faith_pd_boxplot <- ggplot(meta_merged, aes(Tissue, faith_pd)) + 
  geom_boxplot(aes(color = Tissue)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("R_output/tissue_pd.png", tissue_faith_pd_boxplot, height = 3, width = 3)

EcoliTreatment_faith_pd_boxplot <- ggplot(meta_merged, aes(EcoliTreatment, faith_pd)) + ##########not working
  geom_boxplot(aes(color = EcoliTreatment)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("R_output/Ecoli_pd.png", EcoliTreatment_faith_pd_boxplot, height = 3, width = 3)

seedsani_faith_pd_boxplot <- ggplot(meta_merged, aes(SeedSanitization, faith_pd)) + 
  geom_boxplot(aes(color = SeedSanitization)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("R_output/seedsani_pd.png", seedsani_faith_pd_boxplot, height = 3, width = 3)

Age_faith_pd_boxplot <- ggplot(meta_merged, aes(AgeAtHarvest, faith_pd)) + 
  geom_boxplot(aes(color = AgeAtHarvest)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("R_output/age_pd.png", Age_faith_pd_boxplot, height = 3, width = 3)

tissue_evenness_boxplot <- ggplot(meta_merged, aes(Tissue, pielou_evenness)) + 
  geom_boxplot(aes(color = Tissue)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou Evenness", x = "") 
ggsave("R_output/tissue_evenness.png", tissue_evenness_boxplot, height = 3, width = 3)

EcoliTreatment_evenness_boxplot <- ggplot(meta_merged, aes(EcoliTreatment, pielou_evenness)) + ##########not working
  geom_boxplot(aes(color = EcoliTreatment)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou Evenness", x = "") 
ggsave("R_output/Ecoli_evenness.png", EcoliTreatment_evenness_boxplot, height = 3, width = 3)

seedsani_evenness_boxplot <- ggplot(meta_merged, aes(SeedSanitization, pielou_evenness)) + 
  geom_boxplot(aes(color = SeedSanitization)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou Evenness", x = "") 
ggsave("R_output/seedsani_evenness.png", seedsani_evenness_boxplot, height = 3, width = 3)

Age_evenness_boxplot <- ggplot(meta_merged, aes(AgeAtHarvest, pielou_evenness)) + 
  geom_boxplot(aes(color = AgeAtHarvest)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou Evenness", x = "") 
ggsave("R_output/age_evenness.png", Age_evenness_boxplot, height = 3, width = 3)

tissue_obs_feat_boxplot <- ggplot(meta_merged, aes(Tissue, observed_features)) + 
  geom_boxplot(aes(color = Tissue)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("R_output/tissue_obsfeat.png", tissue_obs_feat_boxplot, height = 3, width = 3)

EcoliTreatment_obs_feat_boxplot <- ggplot(meta_merged, aes(EcoliTreatment, observed_features)) + ##########not working
  geom_boxplot(aes(color = EcoliTreatment)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("R_output/Ecoli_obsfeat.png", EcoliTreatment_obs_feat_boxplot, height = 3, width = 3)

seedsani_obs_feat_boxplot <- ggplot(meta_merged, aes(SeedSanitization, observed_features)) + 
  geom_boxplot(aes(color = SeedSanitization)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("R_output/seedsani_obsfeat.png", seedsani_obs_feat_boxplot, height = 3, width = 3)

Age_obs_feat_boxplot <- ggplot(meta_merged, aes(AgeAtHarvest, observed_features)) + 
  geom_boxplot(aes(color = AgeAtHarvest)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("R_output/age_obsfeat.png", Age_obs_feat_boxplot, height = 3, width = 3)

tissue_shannon_boxplot <- ggplot(meta_merged, aes(Tissue, shannon_entropy)) + 
  geom_boxplot(aes(color = Tissue)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon", x = "") 
ggsave("R_output/tissue_shannon.png", tissue_shannon_boxplot, height = 3, width = 3)

EcoliTreatment_shannon_boxplot <- ggplot(meta_merged, aes(EcoliTreatment, shannon_entropy)) + ##########not working
  geom_boxplot(aes(color = EcoliTreatment)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon", x = "") 
ggsave("R_output/Ecoli_shannon.png", EcoliTreatment_shannon_boxplot, height = 3, width = 3)

seedsani_shannon_boxplot <- ggplot(meta_merged, aes(SeedSanitization, shannon_entropy)) + 
  geom_boxplot(aes(color = SeedSanitization)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon", x = "") 
ggsave("R_output/seedsani_shannon.png", seedsani_shannon_boxplot, height = 3, width = 3)

Age_shannon_boxplot <- ggplot(meta_merged, aes(AgeAtHarvest, shannon_entropy)) + 
  geom_boxplot(aes(color = AgeAtHarvest)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon", x = "") 
ggsave("R_output/age_shannon.png", Age_shannon_boxplot, height = 3, width = 3)

##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

##############################################not used in project#############################
# **Normally distributed metrics**

# for continuous variables, we run a general linear model. 
#The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

#glm.evenness.time = glm(pielou_evenness ~ days.since.experiment.start, data=meta)
#summary(glm.evenness.time)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

#plot(pielou_evenness ~ days.since.experiment.start, data=meta)
#Add the glm best fit line
#plot(pielou_evenness ~ days.since.experiment.start, data=meta) + abline(glm.evenness.time)
########################################################################################

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.
gaussian.faith.tissue = glm(faith_pd ~ Tissue, data=meta_merged, family="gaussian")
plot(gaussian.faith.tissue, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.tissue = glm(faith_pd ~ Tissue, data=meta_merged, family="quasipoisson")
plot(qp.faith.tissue, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.tissue)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and days.since.experiment.start.

#Plot
plot(log(faith_pd) ~ Tissue, data=meta_merged, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ Tissue, data=meta_merged, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.tissue)

