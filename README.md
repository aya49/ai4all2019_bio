# ai4all2019_bio

2 week summer enrichment program for grad 10 & 11 girls

links
- drive: https://drive.google.com/drive/folders/1YmWcwDZ4f5B1ajF90gEeKbE60hN1-uBL?usp=sharing
- challenge site: https://www.synapse.org/#!Synapse:syn18380862
- official site: https://www.sfu.ca/computing/inventthefuture.html

## schedule (PST)

meetings
- 2019-05-09 15:00 mentor's startup meeting ([ppt](ITF2019-MentorWelcome.pptx), 
[pdf](ITF2019-MentorWelcome.pdf), [project ideas](https://sfu-db.github.io/bigdata-cmpt733/final-project-sp19.html))
- 2019-05-17 11:00 project pitch ([bio ideas](https://docs.google.com/document/d/1v7Q5Cw732rBZHirZqWpQawUZO749UbMYdlv1ElbI2ZI/edit?usp=sharing))
- 2019-05-22 10:00 pregnancy project webinar https://zoom.us/j/812397944 ([more info](https://www.synapse.org/#!Synapse:syn18380862/discussion/threadId=5365))
- 2019-05-28 15:30 refined project pitch
- 2019-06-28 project code ready
- 2019-07-03 test code on laptops
- 2019-07-10 install environment

program
- 2019-07-19 afternoon project group preferences
- 2019-07-{22:26}\24 project time & event

challenge
- 2019-05-04 sub-challenge 1 open
- 2019-05-22 sub-challenge 1 leaderboard (max 5 submissions)
- 2019-08-15 sub-challenge 1 deadline for submission; 2019-09-13 results announced, 219-11-04 RSG with DREAM conference
- 2019-08-15 sub-challenge 2 open + leaderboard (max 5 submissions)
- 2019-12-05 sub-challenge 2 deadline for submission; 2020-01-03 results announced

data: sub-challenge 1
- input: ```HTA20_RMA/_probeset.RData```
- output: gestational age in weeks for each sample
- ```HTA20_RMA/_probeset.RData```: whole blood preprocessed eset_HTA20/_probeset **(32830/925032 gene/probeset x 367 train + 368 test samples)** expression matrix
  - rownames: ENTREZ-gene(except for “_at” suffix)/probeset IDs
  - colnames: SampleID
- ```anoSC1_v11_nokey.csv```: sample annotation file with the following columns
  - SampleID: unique identifier of the sample (matching the name of the .CEL file in HTA20 folder, except for extension .CEL);
  - GA: gestational age as determined by the last menstrual period and or ultrasound; 
  - Batch: the batch identifier; 
  - Set: name of the source dataset; 
  - Train: 1 for samples to be used for training, 0 for samples to be used for test; 
  - Platform: gene expression platform used to generate the cell files.
- ```preprocess_data_SC1.R```: generates gene level and probeset (exon and exon-junction) level data HTA20_RMA/_probeset.RData from raw (.CEL) files
  - obtained using the RMA implemented in the oligo package of Bioconductor
  - probsets based on a chip definition file generated from the ```pd.hta20.hs.entrezg_0.0.1.tar.gz``` annotation data from http://brainarray.mbni.med.umich.edu; annotation package pd.hta.2.0 (v 3.12.2) from Bioconductor can be used to map probesets to gene identifiers
  - processed separately for each of the 11 experimental batch profiling train + test samples being profiled in each batch
  - then combined across batches 
  - quantile normalized
  - batch effect removal using the removeBatchEffects from the limma package of Bioconductor
- ```pd.hta20.hs.entrezg_0.0.1.tar.gz```: annotation package used to generate gene level expression data.
- ```HT20``` folder (not uploaded 47.36gb): 735 raw .CEL files from the HTA 2.0 platform
