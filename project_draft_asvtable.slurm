#!/bin/sh -l

#SBATCH --nodes=1                   ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 50		         ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name class_project ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

cd /scratch/bell/kquinter/ClassProject

pwd

echo "loading bioinfo"
module load bioinfo

echo "loading qiime2-2023.9"
module load Qiime/2-2023.9

echo "Data Import"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ManifestFile.txt \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path paired-end-demux.qza

echo "making paired-end-demux.qzv"

qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 190 \
  --p-trunc-len-r 180 \
  --p-n-threads 50 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#Remove mitochondria and chloroplast

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

#Rename files so that the filtered one is table.qza

mv table.qza table-with-mitochondria-and-chloroplast.qza
mv table-no-mitochondria-no-chloroplast.qza table.qza

echo "making stats-dada2"

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization stats-dada2.qzv

echo "making table.qzv"

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file MetaData_complete2.txt

echo "making rep-seqs.qzv"

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"
