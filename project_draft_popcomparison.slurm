#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 50                    ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name class_project ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd /scratch/bell/kquinter/ClassProject


pwd

#step 2, Load Modules

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2023.9"
module load Qiime/2-2023.9

#step 3, Run the code you want

wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2024.2/common/gg-13-8-99-515-806-nb-classifier.qza"

echo "Classify and get .qzv"
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


echo "Make a stacked barplot"
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization taxa-bar-plots.qzv


echo "Filtering away some sample types"

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-where "[Tissue]='PlantLeaves'" \
  --o-filtered-table leaves-table.qza

qiime feature-table filter-samples \
  --i-table leaves-table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-where "[AgeAtHarvest]='28 days'" \
  --o-filtered-table 28day-leaves-table.qza

echo "Doing ANCOM"

qiime composition add-pseudocount \
  --i-table 28day-leaves-table.qza \
  --o-composition-table comp-28day-leaves-table.qza

qiime composition ancom \
  --i-table comp-28day-leaves-table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization ancom-SeedSanitization.qzv
  
qiime composition ancom \
  --i-table comp-28day-leaves-table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization ancom-EcoliTreatment.qzv

echo "Collapsing levels"
qiime taxa collapse \
  --i-table 28day-leaves-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table 28day-leaves-table-l5.qza

qiime taxa collapse \
  --i-table 28day-leaves-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table 28day-leaves-table-l6.qza

qiime taxa collapse \
  --i-table 28day-leaves-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table 28day-leaves-table-l7.qza

qiime composition add-pseudocount \
  --i-table 28day-leaves-table-l5.qza \
  --o-composition-table comp-28day-leaves-table-l5.qza

qiime composition add-pseudocount \
  --i-table 28day-leaves-table-l6.qza \
  --o-composition-table comp-28day-leaves-table-l6.qza

qiime composition add-pseudocount \
  --i-table 28day-leaves-table-l7.qza \
  --o-composition-table comp-28day-leaves-table-l7.qza

echo "ANCOM again"

qiime composition ancom \
  --i-table comp-28day-leaves-table-l5.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization l5-ancom-SeedSanitization.qzv

qiime composition ancom \
  --i-table comp-28day-leaves-table-l6.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization l6-ancom-SeedSanitization.qzv

qiime composition ancom \
  --i-table comp-28day-leaves-table-l7.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization l7-ancom-SeedSanitization.qzv

qiime composition ancom \
  --i-table comp-28day-leaves-table-l5.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization l5-ancom-EcoliTreatment.qzv

qiime composition ancom \
  --i-table comp-28day-leaves-table-l6.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization l6-ancom-EcoliTreatment.qzv

qiime composition ancom \
  --i-table comp-28day-leaves-table-l7.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization l7-ancom-EcoliTreatment.qzv

echo "Doing ANCOMBC"

qiime composition ancombc \
  --i-table 28day-leaves-table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'SeedSanitization' \
  --o-differentials ancombc-SeedSanitization.qza
  
qiime composition ancombc \
  --i-table 28day-leaves-table.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'EcoliTreatment' \
  --o-differentials ancombc-EcoliTreatment.qza

qiime composition da-barplot \
  --i-data ancombc-SeedSanitization.qza \
  --p-significance-threshold 0.1 \
  --o-visualization da-barplot-SeedSanitization.qzv
  
qiime composition da-barplot \
  --i-data ancombc-EcoliTreatment.qza \
  --p-significance-threshold 0.1 \
  --o-visualization da-barplot-EcoliTreatment.qzv

echo "Collapsed table made. ANCOMBC again, last step"

qiime composition ancombc \
  --i-table 28day-leaves-table-l5.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'SeedSanitization' \
  --o-differentials l5-ancombc-SeedSanitization.qza

qiime composition ancombc \
  --i-table 28day-leaves-table-l6.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'SeedSanitization' \
  --o-differentials l6-ancombc-SeedSanitization.qza

qiime composition ancombc \
  --i-table 28day-leaves-table-l7.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'SeedSanitization' \
  --o-differentials l7-ancombc-SeedSanitization.qza

qiime composition ancombc \
  --i-table 28day-leaves-table-l5.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'EcoliTreatment' \
  --o-differentials l5-ancombc-EcoliTreatment.qza

qiime composition ancombc \
  --i-table 28day-leaves-table-l6.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'EcoliTreatment' \
  --o-differentials l6-ancombc-EcoliTreatment.qza

qiime composition ancombc \
  --i-table 28day-leaves-table-l7.qza \
  --m-metadata-file MetaData_complete2.txt \
  --p-formula 'EcoliTreatment' \
  --o-differentials l7-ancombc-EcoliTreatment.qza

qiime composition da-barplot \
  --i-data l5-ancombc-SeedSanitization.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l5-da-barplot-SeedSanitization.qzv
 
qiime composition da-barplot \
  --i-data l6-ancombc-SeedSanitization.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot-SeedSanitization.qzv

qiime composition da-barplot \
  --i-data l7-ancombc-SeedSanitization.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l7-da-barplot-SeedSanitization.qzv

qiime composition da-barplot \
  --i-data l5-ancombc-EcoliTreatment.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l5-da-barplot-EcoliTreatment.qzv

qiime composition da-barplot \
  --i-data l6-ancombc-EcoliTreatment.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot-EcoliTreatment.qzv

qiime composition da-barplot \
  --i-data l7-ancombc-EcoliTreatment.qza \
  --p-significance-threshold 0.1 \
  --p-level-delimiter ';' \
  --o-visualization l7-da-barplot-EcoliTreatment.qzv


echo "End time"
date +"%d %B %Y %H:%M:%S"
