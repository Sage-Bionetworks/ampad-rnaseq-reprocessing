###### Rosmap Meta Data Cleaning and Counts combining
# Metadata and counts for sageseqr pipeline are assembled and 
# pushed to synapse
###### Rosmap Meta Data Cleaning and Counts combining
# Metadata and counts for sageseqr pipeline are assembled and 
# pushed to synapse
library(readxl)
library(synapser)
library(data.table)
library(dplyr)
library(tidyr)
library(digest)
library(stringr)


synapser::synLogin()
# Track files for provenance
synids_used <- NULL

##############################################################################
# Combine gene counts files for cohort

# Counts
counts_one <- 'syn22283382'
counts_two <- 'syn22301601'
counts_three <- 'syn22314230'
synids_used <- c(synids_used, counts_one, counts_two, counts_three)

counts_one <- data.table::fread(
  synapser::synGet(counts_one)$path,
  sep = '\t',
  header = T )
counts_two <- data.table::fread(
  synapser::synGet(counts_two)$path,
  sep = '\t',
  header = T )
counts_three <- data.table::fread(
  synapser::synGet(counts_three)$path,
  sep = '\t',
  header = T,
)

counts <- dplyr::left_join(counts_one, counts_two, by = 'feature') %>%
  dplyr::left_join(counts_three, by = 'feature') %>%
  tibble::column_to_rownames(var = "feature")

counts <- counts[ row.names(counts)[ 
  !(row.names(counts) %in% 
      c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous')
  )
],
]
##############################################################################
# Combine MetaData files for cohort

# ID Translation Key:
id_key <- 'syn20818933'
synids_used <- c(synids_used,id_key)
id_key <- read.csv(synapser::synGet(id_key)$path, header = F, stringsAsFactors = F)
id_key$V1 <- stringr::str_pad(as.character(id_key$V1), 8, pad = "0")
id_key <- as.data.frame(id_key[,1:2])
colnames(id_key) <- c('projid','individualID')

# Complete Ages Info (xslx)
ages_a <- 'syn18693152'
synids_used <- c(synids_used, ages_a)
ages_a <- readxl::read_excel(synapser::synGet(ages_a)$path)

# Only one ID from the ID key is missing in ages_a: 
# This indv has no RNA-Seq sample
id_key[ !(id_key$projid %in% ages_a$projid), ]
# projid   patid
# 90214403 R5615105

# Race and Ethnicity Info: 1 = Male, 0 = Female
ages_b <- 'syn18632500'
synids_used <- c(synids_used, ages_b)
ages_b <- readxl::read_excel(synapser::synGet(ages_b)$path) 

# Sequencing Statistics
metrics_one <- 'syn22283384'
metrics_two <- 'syn22301603'
metrics_three <- 'syn22314232'
synids_used <- c(synids_used, metrics_one, metrics_two, metrics_three)

metrics_one <- as.matrix(data.table::fread(
  synapser::synGet(metrics_one)$path,
  sep = '\t',
  header = T))
metrics_two <- as.matrix(data.table::fread(
  synapser::synGet(metrics_two)$path,
  sep = '\t',
  header = T ))
metrics_three <- as.matrix(data.table::fread(
  synapser::synGet(metrics_three)$path,
  sep = '\t',
  header = T ))

metrics <- as.data.frame(rbind(metrics_one, metrics_two, metrics_three), 
                         stringsAsFactors = F
)

for(i in c(2,4:dim(metrics)[2])) {
  metrics[,i] <- as.numeric(metrics[,i])
}

# Clinical Data
clinical <- 'syn3191087'
synids_used <- c(synids_used, clinical)
clinical <- read.csv(synapser::synGet(clinical)$path, stringsAsFactors = F)
clinical$projid <- stringr::str_pad(as.character(clinical$projid), 8, pad = "0")
table( paste0(clinical$projid, '_', clinical$individualID) %in% 
         paste0(id_key$projid, '_', id_key$individualID) 
)

# biospecimin
biospecimin <- 'syn21323366'
ver <- 6
synids_used <- c(synids_used, biospecimin)
biospecimin <- read.csv(synapser::synGet(biospecimin, version=ver)$path, 
                        stringsAsFactors = F
)

biospecimin <- biospecimin[
  biospecimin$assay == 'rnaSeq' & 
    biospecimin$tissue != 'blood',
]

# Assay Data 
assay <- 'syn21088596'
synids_used <- c(synids_used, assay)
assay <- read.csv(synapser::synGet(assay)$path,
                  stringsAsFactors = F
)

#table(assay$specimenID %in% biospecimin$specimenID)
#table(assay[
#    assay$specimenID %in% biospecimin$specimenID, 
#  ]$libraryPreparationMethod
#)
#table(table(assay[
#      assay$specimenID %in% biospecimin$specimenID, 
#    ]$specimenID
#  )
#)

# Annotate everything by sampleID
colnames(metrics)[colnames(metrics) == 'sample'] <- 'specimenID'
table(metrics$specimenID %in% biospecimin$specimenID)

# Combine Seq Metrics and Biospecimin IDs
comb <- dplyr::right_join(biospecimin, metrics, by = 'specimenID')
dim(comb[ is.na(comb$individualID), ])
table(comb[ is.na(comb$individualID), ]$excludeReason)
dim(comb)

# Add in Assay metrics 
comb <- dplyr::right_join(comb, 
                          assay[assay$specimenID %in% comb$specimenID,],
                          by = 'specimenID'
)
dim(comb[ is.na(comb$individualID), ])
table(comb[ is.na(comb$individualID), ]$excludeReason)
dim(comb)

comb <- comb[ is.na(comb$exclude), ]
# 2 weird samples: These samples have no associated individual or patient IDs
# indivdual ID's exist in syn4259690, but these indvidualIDs are missing in 
#   All other files
# Sample IDs: c('800_130701', '764_130520')
synids_used <- c(synids_used, 'syn4259690')
trans_fix <- read.csv(synapser::synGet('syn4259690')$path, stringsAsFactors = F)
trans_fix[ trans_fix$mrna_id %in% c('764_130520','800_130701'),]$projid
# Project IDs: 59204827 and 86956481
c('86956481','59204827') %in% id_key$projid
c('86956481','59204827') %in% clinical$projid

# Removed these 2 samples along with the sample swaps
biospecimin[ biospecimin$specimenID %in% c('800_130701','764_130520'),]
assay[ assay$specimenID %in% c('800_130701','764_130520'),]
clinical[ clinical$specimenID %in% c('800_130701','764_130520'),]

#Combine clinical the sample swaps and nonsense samples are removed
comb <- dplyr::right_join(comb,
                          clinical[ clinical$individualID %in% comb$individualID,], 
                          by = 'individualID')

# Combine and remove 68 samples swaps + 2 samples with indvID issues 
dim(comb[ !is.na(comb$individualID), ])
removed_samples <- colnames(counts)[!(c(colnames(counts)) %in% comb$specimenID)]
comb_sink <- comb

##############################################################################
#Finish Cleaning metadata:
comb$isStranded[ comb$isStranded == 'True'] <- 'TRUE'

remove_na <- NULL
length_issue <- NULL
all_zeros <- NULL
for(name in colnames(comb)){
  tab <- table(is.na(comb[,name]))
  if('TRUE' %in% names(tab)){
    if( as.numeric(table(is.na(comb[,name]))['TRUE']) == dim(comb)[1]) {
      remove_na <- c(remove_na,name)
    }
  }
  tab <- table(comb[,name])
  if('0' %in% names(tab)){
    if(as.numeric(tab['0']) == dim(comb)[1]) {
      all_zeros <- c(all_zeros,name)
    }
  }
  if(length(tab)==1){
    length_issue <- c(length_issue, name)
  }
}

# Variables that are all NAs
message('Variables that are all NAs:')
message(paste0('Variables that are all NAs: ', 
               paste(remove_na,collapse = ', ')
)
)
# specimenIdSource, samplingDate, BrodmannArea, tissueWeight, tissueVolume, 
# fastingState, isPostMortem, samplingAge, visitNumber, exclude, excludeReason, 
# rnaBatch

# Variables that are all Zeros
message(paste0('Variables that are all Zeros: ',
               paste(all_zeros,collapse = ', ')
)
)
# AlignmentSummaryMetrics__BAD_CYCLES, 
# AlignmentSummaryMetrics__PF_HQ_MEDIAN_MISMATCHES, 
# AlignmentSummaryMetrics__PF_NOISE_READS, RnaSeqMetrics__IGNORED_READS

message(paste0('Variables that are all the same value: ',
               paste(length_issue,collapse = ', ')
)
)
# organ, sampleStatus, nucleicAcidSource, cellType, assay, 
# AlignmentSummaryMetrics__BAD_CYCLES, AlignmentSummaryMetrics__CATEGORY, 
# AlignmentSummaryMetrics__PCT_PF_READS, runType, 'isStranded',
# AlignmentSummaryMetrics__PF_HQ_MEDIAN_MISMATCHES, platform, readStrandOrigin,
# AlignmentSummaryMetrics__PF_NOISE_READS, RnaSeqMetrics__IGNORED_READS 

comb <- comb[,
             colnames(comb)[!(colnames(comb) %in% c(remove_na,length_issue,all_zeros))]
]
##############################################################################
# Read Length Issues for Will
# Find the samples that have unexpected mean lengths: sent to Will
##:## Mean_125 <- comb_clean[ comb_clean$AlignmentSummaryMetrics__MEAN_READ_LENGTH == 125, ]$specimenID
##:## Mean_GT_130 <- comb_clean[ comb_clean$AlignmentSummaryMetrics__MEAN_READ_LENGTH > 130, ]$specimenID
##:## table( grepl('RISK_' ,Mean_GT_130) )
##:## table( grepl('RISK_' ,comb_clean$specimenID) )

##:## dc3_prov <- read.csv(synapser::synGet('syn22314227')$path, stringsAsFactors = F)
##:## dc3_flagged <- dc3_prov[ dc3_prov$specimenID %in% Mean_GT_130,]
##:## write.csv(dc3_flagged, file='~/Desktop/dc3_flagged.csv', row.names=F)

##:## dc2_prov <- read.csv(synapser::synGet('syn22301598')$path, stringsAsFactors = F)
##:## dc2_flagged <- dc2_prov[ dc2_prov$specimenID %in% Mean_125,]
##:## write.csv(dc2_flagged, file='~/Desktop/dc2_flagged.csv', row.names=F)
##############################################################################
# Recode variables
table(comb$tissue)

#tissue
regions <- c('DLPFC','ACC','PCC')
names(regions) <- c( 'dorsolateral prefrontal cortex', 
                     'Head of caudate nucleus', 'posterior cingulate cortex' )
comb$tissue <- as.character(regions[comb$tissue])

#msex
table(comb$msex)
sex_trans <- c('female','male')
names(sex_trans) <- c('0','1') 
comb$sex <- sex_trans[as.character(comb$msex)]

comb <- comb[,colnames(comb)[colnames(comb)!='msex']]

#apoe4
apoe_genos <- names(table(comb$apoe_genotype))
apoe_trans <- rep(0, length(apoe_genos))
names(apoe_trans) <- apoe_genos
apoe_trans[ grepl(4, names(apoe_trans)) ] <- 1
apoe_trans[ grepl(44, names(apoe_trans)) ] <- 2
comb$apoe4_allele <- apoe_trans[as.character(comb$apoe_genotype)]

# notes (DataCuts and sequencingBatch -> final batch variable)
#Batch IDs
#Collapse weird "0, 6, 7" into a single number batch
comb[ comb$libraryBatch %in% "0, 6, 7", ]$sequencingBatch <- 9
comb[ comb$libraryBatch %in% "0, 6, 7", ]$libraryBatch <- 9

# Remove letters from data cuts one and two
comb$sequencingBatch <- gsub('NYGC', '', comb$sequencingBatch)
comb$sequencingBatch <- gsub('RISK_', '',comb$sequencingBatch)

colnames(comb)[colnames(comb)=='notes'] <- 'data_contribution'
comb$data_contribution <- gsub('data contribution batch ', '', comb$data_contribution)

comb$final_batch <- paste0(comb$data_contribution, '_', comb$sequencingBatch)

# Diagnosis
comb$diagnosis <- 'OTHER'
comb[ comb$braaksc >= 4 &
        comb$ceradsc <= 2 &
        comb$cogdx == 4  , ]$diagnosis <- 'AD'
comb[ comb$braaksc <= 3 &
        comb$ceradsc >= 3 &
        comb$cogdx == 1  , ]$diagnosis <- 'CT'

colnames(comb) <- gsub('__', '_', colnames(comb))

comb$Tissue.APOE4 <- paste0(comb$tissue,'.', comb$apoe4_allele)
comb$RIN2 <- comb$RIN^2
comb <- comb[, c("individualID", "specimenID", 'projid', 'Study', "tissue",
                 'diagnosis', 'apoe_genotype', 'apoe4_allele', 'Tissue.APOE4', 
                 'pmi', 'braaksc', 'ceradsc', 'cogdx', 'dcfdx_lv', 'sex', 
                 'educ', 'race', 'spanish', 'age_at_visit_max', 'age_first_ad_dx',
                 'age_death', 'cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'RIN', 
                 'RIN2', 'libraryBatch', 'sequencingBatch', 'data_contribution', 
                 'final_batch', 'libraryPrep', 'libraryPreparationMethod', 
                 'readLength', "AlignmentSummaryMetrics_MEAN_READ_LENGTH", 
                 "AlignmentSummaryMetrics_PCT_ADAPTER", 
                 "AlignmentSummaryMetrics_PCT_CHIMERAS",
                 "AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED", 
                 "AlignmentSummaryMetrics_PCT_READS_ALIGNED_IN_PAIRS",
                 "AlignmentSummaryMetrics_PF_ALIGNED_BASES", 
                 "AlignmentSummaryMetrics_PF_HQ_ALIGNED_BASES", 
                 "AlignmentSummaryMetrics_PF_HQ_ALIGNED_Q20_BASES", 
                 "AlignmentSummaryMetrics_PF_HQ_ALIGNED_READS",
                 "AlignmentSummaryMetrics_PF_HQ_ERROR_RATE", 
                 "AlignmentSummaryMetrics_PF_INDEL_RATE",
                 "AlignmentSummaryMetrics_PF_MISMATCH_RATE", 
                 "AlignmentSummaryMetrics_PF_READS",
                 "AlignmentSummaryMetrics_PF_READS_ALIGNED", 
                 "AlignmentSummaryMetrics_READS_ALIGNED_IN_PAIRS",
                 "AlignmentSummaryMetrics_STRAND_BALANCE", 
                 "AlignmentSummaryMetrics_TOTAL_READS",
                 "RnaSeqMetrics_CODING_BASES", 
                 "RnaSeqMetrics_CORRECT_STRAND_READS",
                 "RnaSeqMetrics_INCORRECT_STRAND_READS",
                 "RnaSeqMetrics_INTERGENIC_BASES",
                 "RnaSeqMetrics_INTRONIC_BASES",
                 "RnaSeqMetrics_MEDIAN_3PRIME_BIAS",
                 "RnaSeqMetrics_MEDIAN_5PRIME_BIAS",
                 "RnaSeqMetrics_MEDIAN_5PRIME_TO_3PRIME_BIAS",
                 "RnaSeqMetrics_MEDIAN_CV_COVERAGE",
                 "RnaSeqMetrics_PCT_CODING_BASES",
                 "RnaSeqMetrics_PCT_CORRECT_STRAND_READS",
                 "RnaSeqMetrics_PCT_INTERGENIC_BASES",
                 "RnaSeqMetrics_PCT_INTRONIC_BASES",
                 "RnaSeqMetrics_PCT_MRNA_BASES",
                 "RnaSeqMetrics_PCT_RIBOSOMAL_BASES",
                 "RnaSeqMetrics_PCT_USABLE_BASES",
                 "RnaSeqMetrics_PCT_UTR_BASES",
                 "RnaSeqMetrics_PF_ALIGNED_BASES",
                 "RnaSeqMetrics_PF_BASES",
                 "RnaSeqMetrics_RIBOSOMAL_BASES",
                 "RnaSeqMetrics_UTR_BASES") ]

##############################################################################
# Add PHI Ages over 90 for normalization
ages_a <- as.data.frame(ages_a)
row.names(ages_a) <- ages_a$projid

comb_uncensored <- comb
comb_uncensored$age_first_ad_dx <-  ages_a[ comb_uncensored$projid,]$age_first_ad_dx
comb_uncensored$age_death <-  ages_a[ comb_uncensored$projid,]$age_death
comb_uncensored$age_at_visit_max <-  ages_a[ comb_uncensored$projid,]$age_at_visit_max
comb_uncensored$cts_mmse30_lv <-  ages_a[ comb_uncensored$projid,]$cts_mmse30_lv
comb_uncensored$cts_mmse30_first_ad_dx <-  ages_a[ comb_uncensored$projid,]$cts_mmse30_first_ad_dx

comb_uncensored <- comb_uncensored[!is.na(comb_uncensored$apoe4_allele),]
comb <- comb[!is.na(comb$apoe4_allele),]
comb_uncensored <- comb_uncensored[!is.na(comb_uncensored$pmi),]
comb <- comb[!is.na(comb$pmi),]
comb_uncensored <- comb_uncensored[!is.na(comb_uncensored$RIN),]
comb <- comb[!is.na(comb$RIN),]

counts <- counts[,comb_uncensored$specimenID]
counts <- counts[,comb_uncensored$specimenID]

#Toss out pipeline isssues
exclude_samples_sex_swap <- c("327_120501", "Sample_R4361022-DLPFC", 
  "Sample_R2730285-DLPFC", "Sample_R8261694-DLPFC", "RISK_333")
exclude_samples_outliers <- c('380_120503', 'RISK_197', 'RISK_278', 'RISK_282', 
                                'RISK_283', 'RISK_295', 'RISK_300', 'RISK_303', 
                                'RISK_305', 'RISK_310', 'RISK_314', 'RISK_63', 
                                'RISK_121_redo', 'RISK_125_redo', 'RISK_171_redo',
                                'RISK_218', 'RISK_243', 'RISK_254', 'RISK_260', 
                                'RISK_269', 'RISK_280', 'RISK_418', 'RISK_81',
                                'RISK_93', 'RISK_321', 'RISK_15', 'RISK_167', 
                              'RISK_29', 'RISK_299', 'RISK_324', 'RISK_329', 
                              'RISK_33', 'RISK_37', 'RISK_439', 'RISK_455', 
                              'RISK_55', 'RISK_59', 'RISK_65', 'RISK_75', 
                              'RISK_79', 'RISK_87', 'RISK_89', 'RISK_101',
                              'RISK_131_redo', 'RISK_137_redo', 'RISK_208',
                              'RISK_217', 'RISK_259', 'RISK_47', 'RISK_97')
  comb_uncensored <- comb_uncensored[ 
    !(comb_uncensored$specimenID %in% 
        c(exclude_samples_sex_swap,exclude_samples_outliers)) ,
  ]
  
  counts <- counts[, 
                   colnames(counts)[
                     !(colnames(counts) %in% 
                       c(exclude_samples_sex_swap,exclude_samples_outliers)
                     )]
                   ]
##############################################################################
#Upload full file and SageSeqr input version to synapse - internal Sage Location:
internal_parentid <- 'syn25784097'
folder_loc <- 'syn25792226'

#Set Activity
#activity <- synapser::synGetEntity('syn25792226')
activity <- synapser::synGet('syn25792226')

##Set Annotations:
all.annotations = list(
  dataType = c('clinical','geneExpression'),
  resourceType = 'metadata',
  metadataType = 'analytical covariates',
  isModelSystem = 'FALSE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'csv',
  grant = 'U01AG046152',
  species = 'Human',
  organ = 'brain',
  tissue = c('dorsolateral prefrontal cortex', 
             'Head of caudate nucleus', 
             'posterior cingulate cortex'
  ),
  study = c('ROSMAP','rnaSeqReprocessing'), 
  consortium = 'AMP-AD',
  assay = 'rnaSeq'
)

#Github Code Pull
thisRepo <- githubr::getRepo(
  repository = "Sage-Bionetworks/ampad-rnaseq-reprocessing",
  ref="branch",
  refName='main'
)
thisFile <- githubr::getPermlink(
  repository = thisRepo,
  repositoryPath=paste0('code/metadata_preprocessing/',
                        'rosmap_preprocessing.R'
  )
)

activityName = 'Full RosMap Metadata'
activityDescription = 'Cleaned Codified and Recoded Ages Uncensored Metadata'

write.csv(comb_uncensored,
          file = 'Full_ROSMAP_RNASeq_Covariates_Uncensored.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Full_ROSMAP_RNASeq_Covariates_Uncensored.csv',
  name = 'RosMap Ages Uncensored Full Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Full_ROSMAP_RNASeq_Covariates_Uncensored.csv")

## Upload Sageseqr file and SageSeqr input version to synapse - internal Sage Location:
# Full Rosmap Cohort:
comb_uncensored_sageseqr <- comb_uncensored[,c(
  'specimenID',	'diagnosis',	'tissue',	'race',	'spanish',	'apoe4_allele',
  'sex',	'final_batch',	'pmi',	'RIN',	'RIN2',	'age_death',	
  'AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED',	
  'RnaSeqMetrics_PCT_INTRONIC_BASES', 'RnaSeqMetrics_PCT_INTERGENIC_BASES',
  'RnaSeqMetrics_PCT_CODING_BASES'
)]

activityDescription = 'Cleaned Codified and Recoded Ages Uncensored Metadata'

write.csv(comb_uncensored_sageseqr,
          file = 'Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv',
  name = 'RosMap Ages Uncensored Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv")

# Indv Tissue Metadata for sage-seqr
for (tiss in c('ACC', 'DLPFC', 'PCC' )) {
  tiss_file <- comb_uncensored_sageseqr[comb_uncensored_sageseqr$tissue == tiss,]
  tiss_file <- tiss_file[,colnames(tiss_file)[!(colnames(tiss_file)%in%'tissue')]]
 
  activityDescription = 'Cleaned Codified and Recoded Ages Uncensored Metadata'
  
  write.csv(tiss_file,
            file = paste0(tiss,
                          '_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'
                        ),
            row.names = F,
            quote = F
  )
  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path= paste0(tiss,'_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'),
    name = paste0(tiss, ' RosMap Ages Uncensored Sageseqr Input Covariates'),
    parentId='syn26242659' ),
    used = synids_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0(tiss,'_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'))
}


# Indv Tissue Metadata for sage-seqr
for (tiss in c('ACC', 'DLPFC', 'PCC' )) {
  tiss_file <- comb_uncensored_sageseqr[comb_uncensored_sageseqr$tissue == tiss,]
  tiss_file <- tiss_file[,colnames(tiss_file)[!(colnames(tiss_file)%in%'tissue')]]
  tiss_file <- tiss_file[ tiss_file$RIN >= 5, ]
  activityDescription = 'Cleaned Codified and Recoded Ages Uncensored Metadata'
  
  write.csv(tiss_file,
            file = paste0(tiss,
                          'Rinover5_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'
            ),
            row.names = F,
            quote = F
  )
  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path= paste0(tiss,'Rinover5_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'),
    name = paste0(tiss, 'Rin over 5 RosMap Ages Uncensored Sageseqr Input Covariates'),
    parentId='syn26242659' ),
    used = synids_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0(tiss,'Rinover5_Sageseqr_ROSMAP_RNASeq_Covariates_Uncensored.csv'))
}

## Upload Sageseqr ages Censored file to synapse
comb_censored_sageseqr <- comb[,colnames(comb_uncensored_sageseqr)]

children <- synapser::synGetChildren(activity$properties$id)$asList()
cesoredmeta_used <- NULL
for(i in 1:length(children)) {
  if(children[[i]]$name == "RosMap Ages Uncensored Sageseqr Input Covariates") {
    cesoredmeta_used <- children[[i]]$id
  }
}

censordmeta_parentid <- 'syn25808143'
activity <- synapser::synGetEntity(censordmeta_parentid)
activityName = 'Sageseqr Input Metadata'
activityDescription = 'Cleaned Codified and Recoded Ages Censored Metadata'

write.csv(comb_censored_sageseqr,
          file = 'Sageseqr_ROSMAP_RNASeq_Covariates_Censored.csv',
          row.names = F,
          quote = F
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='Sageseqr_ROSMAP_RNASeq_Covariates_Censored.csv',
  name = 'RosMap Ages Censored Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = cesoredmeta_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_ROSMAP_RNASeq_Covariates_Censored.csv")

## Upload Sageseqr counts
counts_used <- c('syn22283382', 'syn22301601', 'syn22314230')
counts_parentid <- 'syn25808173'
activity <- synapser::synGet(counts_parentid)
activityName = 'Sageseqr Input Counts'
activityDescription = 'Combined ROSMAP RNASeq Counts for input to SageSeqr'

all.annotations.expression = list(
  dataType = 'geneExpression',
  dataSubtype = 'processed',
  resourceType = 'experimentalData',
  assay = 'rnaSeq',
  nucleicAcidSource = "bulk cell",
  isModelSystem = 'FALSE',
  isConsortiumAnalysis = 'TRUE',
  isMultiSpecimen = 'TRUE',
  fileFormat = 'txt',
  grant = 'U01AG046152',
  species = 'Human',
  organ = 'brain',
  tissue = c('dorsolateral prefrontal cortex', 
             'head of caudate nucleus', 
             'posterior cingulate cortex'
  ),
  study = c('ROSMAP','rnaSeqReprocessing'), 
  consortium = 'AMP-AD'
)

counts_write <- counts
counts_write <- counts_write[ ,comb_uncensored_sageseqr$specimenID ]
counts_write$feature <- row.names(counts_write)
counts_write <- counts_write[,c('feature',
                                colnames(counts_write)[!(colnames(counts_write) %in% 'feature')]
)
]


write.table(counts_write,
            file = 'ROSMAP_counts.txt',
            row.names = F,
            col.names = T,
            quote = F,
            sep = '\t'
)

ENRICH_OBJ <- synapser::synStore( synapser::File( 
  path='ROSMAP_counts.txt',
  name = 'RosMap Sageseqr Input Counts',
  parentId=activity$properties$id ),
  used = counts_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
file.remove("ROSMAP_counts.txt")

anno_trans <- c('dorsolateral prefrontal cortex', 
  'head of caudate nucleus', 
  'posterior cingulate cortex'
)
names(anno_trans) <- c('DLPFC','ACC', 'PCC' )

for (tiss in c('ACC', 'DLPFC', 'PCC' )) {
  tiss_file <- comb_uncensored_sageseqr[comb_uncensored_sageseqr$tissue == tiss,]
  tiss_file <- tiss_file[,colnames(tiss_file)[!(colnames(tiss_file)%in%'tissue')]]
  
  
  counts_used <- c('syn22283382', 'syn22301601', 'syn22314230')
  counts_parentid <- 'syn25808173'
  activity <- synapser::synGet(counts_parentid)
  activityName = paste0(tiss, ' Sageseqr Input Counts')
  activityDescription = paste0( tiss, ' ROSMAP RNASeq Counts for input to SageSeqr')
  
  all.annotations.expression = list(
    dataType = 'geneExpression',
    dataSubtype = 'processed',
    resourceType = 'experimentalData',
    assay = 'rnaSeq',
    nucleicAcidSource = "bulk cell",
    isModelSystem = 'FALSE',
    isConsortiumAnalysis = 'TRUE',
    isMultiSpecimen = 'TRUE',
    fileFormat = 'txt',
    grant = 'U01AG046152',
    species = 'Human',
    organ = 'brain',
    tissue = as.character(anno_trans[tiss]),
    study = c('ROSMAP','rnaSeqReprocessing'), 
    consortium = 'AMP-AD'
  )
  
  counts_write <- counts
  counts_write <- counts_write[ ,comb_uncensored_sageseqr$specimenID ]
  counts_write$feature <- row.names(counts_write)
  
  counts_write <- counts_write[,c('feature',
                                  colnames(counts_write)[colnames(counts_write) %in% tiss_file$specimenID]
  )
  ]
  
  
  write.table(counts_write,
              file = paste0(tiss, '_ROSMAP_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path=paste0(tiss, '_ROSMAP_counts.txt'),
    name = paste0(tiss,' RosMap Sageseqr Input Counts'),
    parentId='syn26242660'),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0(tiss, '_ROSMAP_counts.txt'))
}



for (tiss in c('ACC', 'DLPFC', 'PCC' )) {
  tiss_file <- comb_uncensored_sageseqr[comb_uncensored_sageseqr$tissue == tiss & 
                                          comb_uncensored_sageseqr$RIN >= 5,]
  tiss_file <- tiss_file[,colnames(tiss_file)[!(colnames(tiss_file)%in%'tissue')]]
  
  
  counts_used <- c('syn22283382', 'syn22301601', 'syn22314230')
  counts_parentid <- 'syn25808173'
  activity <- synapser::synGet(counts_parentid)
  activityName = paste0( tiss, ' Rin over 5 Sageseqr Input Counts')
  activityDescription = paste0( tiss, ' Rin over 5 Combined ROSMAP RNASeq Counts for input to SageSeqr')
  
  all.annotations.expression = list(
    dataType = 'geneExpression',
    dataSubtype = 'processed',
    resourceType = 'experimentalData',
    assay = 'rnaSeq',
    nucleicAcidSource = "bulk cell",
    isModelSystem = 'FALSE',
    isConsortiumAnalysis = 'TRUE',
    isMultiSpecimen = 'TRUE',
    fileFormat = 'txt',
    grant = 'U01AG046152',
    species = 'Human',
    organ = 'brain',
    tissue = as.character(anno_trans[tiss]),
    study = c('ROSMAP','rnaSeqReprocessing'), 
    consortium = 'AMP-AD'
  )
  
  counts_write <- counts
  counts_write <- counts_write[ ,comb_uncensored_sageseqr$specimenID ]
  counts_write$feature <- row.names(counts_write)
  
  counts_write <- counts_write[,c('feature',
                                  colnames(counts_write)[colnames(counts_write) %in% tiss_file$specimenID]
  )
  ]
  
  
  write.table(counts_write,
              file = paste0(tiss, '_Rinover5_ROSMAP_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File( 
    path=paste0(tiss, '_Rinover5_ROSMAP_counts.txt'),
    name = paste0(tiss,' Rin over 5 RosMap Sageseqr Input Counts'),
    parentId='syn26242660'),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0(tiss, '_Rinover5_ROSMAP_counts.txt'))
}
 






