#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 50		         ##number of requested CPUs
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

echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

echo "Rarefaction"

#consult table.qzv for max-depth
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1300 \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization alpha-rarefaction.qzv

echo "Calculate diversity metrics." 
#You will need to designate a sampling depth. You might want to consult the rarifaction output prior to this step.

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1250 \ 
  --m-metadata-file MetaData_complete2.txt \
  --output-dir core-metrics-results

#Statistical testing
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column Tissue \
  --o-visualization core-metrics-results/unweighted-unifrac-tissue-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization core-metrics-results/unweighted-unifrac-treatment-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization core-metrics-results/unweighted-unifrac-sanitization-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column AgeAtHarvest \
  --o-visualization core-metrics-results/unweighted-unifrac-age-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column Tissue \
  --o-visualization core-metrics-results/weighted-unifrac-tissue-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column EcoliTreatment \
  --o-visualization core-metrics-results/weighted-unifrac-treatment-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column SeedSanitization \
  --o-visualization core-metrics-results/weighted-unifrac-sanitization-significance.qzv \
  --p-pairwise
  
  qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file MetaData_complete2.txt \
  --m-metadata-column AgeAtHarvest \
  --o-visualization core-metrics-results/weighted-unifrac-age-significance.qzv \
  --p-pairwise
  
  #Make ordination
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-AgeAtHarvest.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file MetaData_complete2.txt \
  --o-visualization core-metrics-results/bray-curtis-emperor-AgeAtHarvest.qzv
  
 echo "End time"
date +"%d %B %Y %H:%M:%S"
