###### Mayo Meta Data Cleaning and Counts combining
# Metadata and counts for sageseqr pipeline are assembled and
# pushed to synapse
library(readxl)
library(synapser)
library(data.table)
library(dplyr)
library(tidyr)
library(digest)
library(stringr)

# Track files for provenance
synids_used <- NULL

### Data download

# Get reprocessed sample counts for temporal cortex
#counst are here: 'syn20825471'
count <- 'syn21544635';
synids_used <- count
count <- synapser::synGet(count)$path %>%
  read.table(header=T, sep='\t', check.names = F, row.names = 1)

count <- count[ row.names(count)[
  !(row.names(count) %in%
      c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous')
  )
],]

# Get clinical metadata
clinical_md <- 'syn23277389'
synids_used = c(synids_used, clinical_md)
clinical_md <- synGet(clinical_md)$path %>%
  read.csv(header=T, stringsAsFactors=F)

#Braak, CERAD, Apoe Geno Standard names
colnames(clinical_md)[ colnames(clinical_md) == "Braak" ] <- 'braaksc'
colnames(clinical_md)[ colnames(clinical_md) == "CERAD" ] <- 'ceradsc'
colnames(clinical_md)[ colnames(clinical_md) == "apoeGenotype" ] <- 'apoe_genotype'
colnames(clinical_md)[ colnames(clinical_md) == "ageDeath" ] <- 'age_death'

# APOE Alleles
clinical_md$apoe4_allele <- NA
clinical_md[clinical_md$apoe_genotype %in% c('33','23','22'),]$apoe4_allele <- 0
clinical_md[clinical_md$apoe_genotype %in% c('34','24'),]$apoe4_allele <- 1
clinical_md[clinical_md$apoe_genotype %in% c('44'),]$apoe4_allele <- 2

#Recoad ethnicity to match  ROSMAP:
#1=White
#2=Black, Negro, African-American
#3=Native American, Indian
#4=Eskimo
#5=Aleut
#6=Asian or Pacific Island
#8=REFUSAL
#9=DON'T KNOW
clinical_md[ is.na(clinical_md$race),]$race <- 9
clinical_md[ clinical_md$race %in% 'White',]$race <- 1

#Code Diagnosis According to re-seq paper: 'AD', 'CT', 'OTHER'
clinical_md[clinical_md$diagnosis %in% 'Alzheimer Disease',]$diagnosis <- 'AD'
clinical_md[clinical_md$diagnosis %in% 'control',]$diagnosis <- 'CT'
clinical_md[
  clinical_md$diagnosis %in% 'progressive supranuclear palsy',
]$diagnosis <- 'PSP'
clinical_md[
  clinical_md$diagnosis %in% 'pathological aging',
]$diagnosis <- 'PATH_AGE'

# Get picard metrics from synapse
metrics <- 'syn21544637';
synids_used = c(synids_used, metrics)
metrics <- synapser::synGet(metrics)$path %>%
  read.table(sep='\t',header=T)

colnames(metrics) <- gsub('__', '_', colnames(metrics))
colnames(metrics)[colnames(metrics) == 'sample'] <- 'specimenID'

## Assay MetaData
assay <- 'syn20827193'
synids_used = c(synids_used, assay)
assay <- synapser::synGet(assay)$path %>%
  read.csv(header=T, stringsAsFactors = F)

## BioSpecimin Metadata
biospecimin <- 'syn20827192'
synids_used = c(synids_used, biospecimin)
biospecimin <- synapser::synGet(biospecimin)$path %>%
  read.csv(header=T, stringsAsFactors = F)
biospecimin <- biospecimin[biospecimin$assay == 'rnaSeq',]

biospecimin[biospecimin$tissue == 'cerebellum',]$tissue <- 'CBE'
biospecimin[biospecimin$tissue == 'temporal cortex',]$tissue <- 'TCX'


# Merge all metadata
metadata = biospecimin %>%
  dplyr::left_join( clinical_md, by = 'individualID' ) %>%
  dplyr::inner_join(assay, by='specimenID') %>%
  dplyr::inner_join(metrics, by='specimenID')

# Samples excluded by MAYO investigators
samples_exclude = c('syn6126119', 'syn6126114')
synids_used = c(synids_used, 'syn6126119', 'syn6126114')
samples_exclude <- synapser::synGet(samples_exclude[1])$path %>%
  read.table(header=T, sep='\t') %>%
  rbind( read.table(
    synapser::synGet(samples_exclude[2])$path,
    header=T,
    sep='\t'
  ))

colnames(count) <- gsub('_CER', '_CBE', colnames(count))
metadata$specimenID <- gsub('_CER', '_CBE', metadata$specimenID)
samples_exclude$Sample.Name <- gsub('_CER', '_CBE', samples_exclude$Sample.Name)

metadata = metadata %>%
  dplyr::filter(specimenID %in% colnames(count)) %>% # No Samples
  dplyr::filter(!is.na(RIN)) %>% # 6 samples (4 missing age_death, all missing source)
  dplyr::filter(!is.na(RnaSeqMetrics_PCT_INTRONIC_BASES)) %>% # No Samples
  dplyr::filter(!is.na(age_death)) %>% # 4 Samples
  dplyr::filter(!is.na(specimenIdSource)) %>% # 6 samples
  dplyr::filter(!(specimenID %in% samples_exclude$Sample.Name)) %>% # 48 samples
  as.data.frame()

# Get RIN square
metadata = metadata %>%
  dplyr::mutate(RIN2 = RIN^2)

metadata$Tissue.APOE4 = paste(metadata$tissue, metadata$apoe4_allele, sep = '.')

# Match covariates to expression data
indToRetain = intersect(metadata$specimenID, colnames(count))
indRemoved = setdiff(colnames(count), metadata$specimenID)
count <- count[, indToRetain]
rownames(metadata) = metadata$specimenID
metadata = metadata[indToRetain,]

#Get gene specific parameters from synapse
gene_param  <- 'syn23625800'
synids_used = c(synids_used, gene_param)
gene_param = synapser::synGet('syn23625800')$path %>%
  read.csv(header=T)

gene_len = dplyr::select(gene_param, ensembl_gene_id, gene_length) %>% unique()
rownames(gene_len) = gene_len$ensembl_gene_id
gene_gc = dplyr::select(gene_param, ensembl_gene_id, percentage_gene_gc_content) %>% unique()
rownames(gene_gc) = gene_gc$ensembl_gene_id

################################################################################
### Predict missing PMI
metadata$individualID <- as.character(metadata$individualID)

# Get PMI
pmi = metadata[,c('pmi', 'specimenID', 'individualID'), drop = F]
pmi$pmi = pmi$pmi

# Get expression (convert counts to cpm)
expr = limma::voom(count, design = NULL)$E
ind = which(rowSums(expr>=1)/dim(expr)[2] >= 0.5)
expr = expr[ind,]

# Combine expression from different regions of brain of an individual into one
expr = plyr::dlply(pmi, ('individualID'), .fun = function(x, expr){
  data.frame(value = rowMeans(expr[,x$specimenID,drop = F], na.rm = T)) %>%
    plyr::rename(c('value' = unique(x$individualID)))
}, expr, .parallel = F) %>%
  do.call(cbind,.)

# Split data into training and prediction set
pmi = dplyr::select(pmi, individualID, pmi) %>% unique()
rownames(pmi) = pmi$individualID
pmi = pmi[colnames(expr),]

# Fit a logistic regression model
x.train = expr[,!is.na(pmi$pmi)]
x.predict = expr[,is.na(pmi$pmi)]
y.train = pmi[!is.na(pmi$pmi),]
y.train = y.train[colnames(x.train), 'pmi']
plr.model = glmnet::cv.glmnet(t(x.train), y.train, family = 'gaussian')
coeff = as.matrix(glmnet::coef.glmnet(plr.model, s = 'lambda.min'))
ssres = sum((predict(plr.model, t(x.train), s = 'lambda.min', type = 'response') - y.train)^2)
sstot = sum((y.train-mean(y.train, na.rm = T))^2, na.rm = T)
R2.allgenes = 1-(ssres/sstot)
y.predict.all = predict(plr.model, t(x.predict), s = 'lambda.min', type = 'response') %>%
  CovariateAnalysis::rownameToFirstColumn('individualID') %>%
  plyr::rename(c('lambda.min' = 'pmi'))

# Get genes that are predictors of pmi and use them as control probes
controlprobes_id = 'syn6145639'
synids_used = c(synids_used, controlprobes_id)
controlprobes_id = fread(synGet(controlprobes_id)$path, data.table = F) %>%
  dplyr::mutate(specimenID = x) %>%
  tidyr::separate(specimenID, c('ensembl_gene_id', 'position'), sep = '\\.')
row.names(expr) <- do.call(rbind, strsplit(rownames(expr), '[.]'))[,1]

# Initial normalisation of gene expression
ind = rownames(expr) %in% controlprobes_id$ensembl_gene_id
expr = expr[ind,]
# Fit a logistic regression model
x.train = expr[,!is.na(pmi$pmi)]
x.predict = expr[,is.na(pmi$pmi)]
y.train = pmi[!is.na(pmi$pmi),]
y.train = y.train[colnames(x.train), 'pmi']
plr.model = glmnet::cv.glmnet(t(x.train), y.train, family = 'gaussian')
coeff = as.matrix(glmnet::coef.glmnet(plr.model, s = 'lambda.min'))
ssres = sum((predict(plr.model, t(x.train), s = 'lambda.min', type = 'response') - y.train)^2)
sstot = sum((y.train-mean(y.train, na.rm = T))^2, na.rm = T)
R2.smallgenes = 1-(ssres/sstot)
y.predict = predict(plr.model, t(x.predict), s = 'lambda.min', type = 'response') %>%
  CovariateAnalysis::rownameToFirstColumn('individualID') %>%
  plyr::rename(c('lambda.min' = 'pmi'))
pmi = pmi %>%
  dplyr::select(individualID, pmi) %>%
  dplyr::filter(!is.na(pmi)) %>%
  list(y.predict) %>%
  data.table::rbindlist(use.names = T, fill = T)
pmi$pmi[pmi$pmi < 0] = 0
metadata = metadata %>%
  dplyr::select(-pmi) %>%
  dplyr::left_join(pmi)
################################################################################
# Remove the factors that don't matter:
remove_na <- NULL
length_issue <- NULL
all_zeros <- NULL
for(name in colnames(metadata)){
  tab <- table(is.na(metadata[,name]))
  if('TRUE' %in% names(tab)){
    if( as.numeric(table(is.na(metadata[,name]))['TRUE']) == dim(metadata)[1]) {
      remove_na <- c(remove_na,name)
    }
  }
  tab <- table(metadata[,name])
  if('0' %in% names(tab)){
    if(as.numeric(tab['0']) == dim(metadata)[1]) {
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

# Variables that are all Zeros
message(paste0('Variables that are all Zeros: ',
               paste(all_zeros,collapse = ', ')
)
)

message(paste0('Variables that are all the same value: ',
               paste(length_issue,collapse = ', ')
)
)

total_metaata <- metadata
metadata <- metadata[,
             colnames(metadata)[!(colnames(metadata) %in% c(remove_na,length_issue,all_zeros))]
]

ordered_columns <- c("individualID", "specimenID", "specimenIdSource", "tissue",
                     "exclude", "sex", "age_death",  "pmi", "diagnosis",
                     "braaksc", "thal", "apoe_genotype", "apoe4_allele",
                     "Tissue.APOE4", "RIN", "RIN2", "flowcell",
                     colnames(metadata)[
                       grepl('AlignmentSummaryMetrics_',colnames(metadata))
                     ],
                     colnames(metadata)[
                       grepl('RnaSeqMetrics_',colnames(metadata))
                     ]
)

metadata <- metadata[,ordered_columns]
total_metaata <- total_metaata[,c(
    ordered_columns,
    colnames(total_metaata)[!(colnames(total_metaata) %in% ordered_columns)]
  )
]
# flowcell
# Un-mask Age Death for sageseqr input

#Remove Outliers:
indToRemove = c("11311_CER", "1923_CER", "1923_TCX", "11396_TCX", "11294_TCX", "11408_TCX",
                '1950_TCX', '1950_CER', '1925_TCX', '1957_CER')

metadata_sageqr = metadata %>%
  dplyr::mutate(age_death = gsub("[+]", "", age_death))

metadata_sageqr <- metadata_sageqr[
  !(metadata_sageqr$specimenID %in% indToRemove),
]

count <- count[, !(colnames(count) %in% indToRemove)]
metadata_sageqr <- metadata_sageqr[ ,
  c( 'specimenID', 'tissue', 'diagnosis',
     'apoe4_allele',	'sex', 'flowcell',	'pmi', 'RIN',
     'RIN2',	'age_death', 'AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED',
     'RnaSeqMetrics_PCT_INTRONIC_BASES', 'RnaSeqMetrics_PCT_INTERGENIC_BASES',
     'RnaSeqMetrics_PCT_CODING_BASES'
    )
]

################################################################################
                    ## --   Push to Synapse   --  ##
#Upload full file and SageSeqr input version to synapse - internal Sage Location:
internal_parentid <- 'syn25832833'
folder_loc <- 'syn25832834'

#Set Activity
activity <- synapser::synGet(folder_loc)

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
  tissue = c('temporal cortex',
             'cerebellum'
  ),
  study = c('MayoRNAseq','rnaSeqReprocessing'),
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
                        'mayo_preprocessing.R'
  )
)

activityName = 'Full Mayo Metadata'
activityDescription = 'Codified and Recoded Metadata'

write.csv(total_metaata,
          file = 'Full_Mayo_RNASeq_Covariates.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='Full_Mayo_RNASeq_Covariates.csv',
  name = 'Mayo Full Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Full_Mayo_RNASeq_Covariates.csv")

## Upload Sageseqr file and SageSeqr input version to synapse - internal Sage Location:
activityDescription = 'Cleaned Codified and Recoded Metadata'

write.csv(metadata,
          file = 'Cleaned_Mayo_RNASeq_Covariates.csv',
          row.names = F,
          quote = F
)
ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='Cleaned_Mayo_RNASeq_Covariates.csv',
  name = 'Mayo Cleaned Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Cleaned_Mayo_RNASeq_Covariates.csv")

## Upload Sageseqr file to synapse
sageseqrmeta_parentid <- 'syn25808143'
activity <- synapser::synGet(sageseqrmeta_parentid)
activityName = 'Sageseqr Input Metadata'
activityDescription = 'Cleaned Codified and Recoded Metadata'

write.csv(metadata_sageqr,
          file = 'Sageseqr_Mayo_RNASeq_Covariates_Censored.csv',
          row.names = F,
          quote = F
)

ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='Sageseqr_Mayo_RNASeq_Covariates_Censored.csv',
  name = 'Mayo Sageseqr Input Covariates',
  parentId=activity$properties$id ),
  used = synids_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
file.remove("Sageseqr_Mayo_RNASeq_Covariates_Censored.csv")

# write the independent tissue metadata files
row.names(metadata) <- metadata$specimenID
for( tis in c('TCX','CBE')){
  meta_write <- metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]
  write.csv(meta_write[ , colnames(meta_write)[!(colnames(meta_write) %in% 'tissue')]],
            file = paste0('Sageseqr_Mayo_', tis,'_RNASeq_Covariates_Censored.csv'),
            row.names = F,
            quote = F
  )

  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('Sageseqr_Mayo_', tis,'_RNASeq_Covariates_Censored.csv'),
    name = paste0('Mayo ', tis,' Sageseqr Input Covariates'),
    parentId=activity$properties$id ),
    used = synids_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0("Sageseqr_Mayo_", tis,"_RNASeq_Covariates_Censored.csv"))
  
  #Braak
  meta_write <- metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]
  meta_write$braaksc <- metadata[meta_write$specimenID,]$braaksc
  meta_write <- meta_write[complete.cases(meta_write),]
  write.csv(meta_write[ , colnames(meta_write)[!(colnames(meta_write) %in% 'tissue')]],
            file = paste0('Sageseqr_Mayo_', tis,'_RNASeq_Braak.csv'),
            row.names = F,
            quote = F
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('Sageseqr_Mayo_', tis,'_RNASeq_Braak.csv'),
    name = paste0('Mayo ', tis,' Sageseqr Input Braak Covariates'),
    parentId=activity$properties$id ),
    used = synids_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0("Sageseqr_Mayo_", tis,"_RNASeq_Braak.csv"))
  
  #Thal
  meta_write <- metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]
  meta_write$thal <- metadata[meta_write$specimenID,]$thal
  meta_write <- meta_write[complete.cases(meta_write),]
  write.csv(meta_write[ , colnames(meta_write)[!(colnames(meta_write) %in% 'tissue')]],
            file = paste0('Sageseqr_Mayo_', tis,'_RNASeq_Thal.csv'),
            row.names = F,
            quote = F
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('Sageseqr_Mayo_', tis,'_RNASeq_Thal.csv'),
    name = paste0('Mayo ', tis,' Sageseqr Input Thal Covariates'),
    parentId=activity$properties$id ),
    used = synids_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations)
  file.remove(paste0("Sageseqr_Mayo_", tis,"_RNASeq_Thal.csv"))
  
}
## Upload Sageseqr counts
counts_used <- c('syn21544635')
counts_parentid <- 'syn25808173'
activity <- synapser::synGet(counts_parentid)
activityName = 'Sageseqr Input Counts'
activityDescription = 'Mayo RNASeq Counts for input to SageSeqr'

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
  tissue = c('temporal cortex',
             'cerebellum'
  ),
  study = c('MayoRNAseq','rnaSeqReprocessing'),
  consortium = 'AMP-AD'
)

counts_write <- count
counts_write$feature <- row.names(counts_write)
counts_write <- counts_write[,c('feature',
                                colnames(counts_write)[!(colnames(counts_write) %in% 'feature')]
)
]
write.table(counts_write,
            file = 'MAYO_counts.txt',
            row.names = F,
            col.names = T,
            quote = F,
            sep = '\t'
)

ENRICH_OBJ <- synapser::synStore( synapser::File(
  path='MAYO_counts.txt',
  name = 'MAYO Sageseqr Input Counts',
  parentId=activity$properties$id ),
  used = counts_used,
  activityName = activityName,
  executed = thisFile,
  activityDescription = activityDescription
)
synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
file.remove("MAYO_counts.txt")


for( tis in c('TCX','CBE')){
  
  keeps <- c('feature', as.character(
    metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]$specimenID
  ))

  write.table(counts_write[,keeps],
              file = paste0('MAYO_', tis,'_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )

  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('MAYO_', tis,'_counts.txt'),
    name = paste0('MAYO ', tis, ' Sageseqr Input Counts'),
    parentId=activity$properties$id ),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0('MAYO_', tis,'_counts.txt'))
  
  #Braaksc
  meta_write <- metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]
  meta_write$braaksc <- metadata[meta_write$specimenID,]$braaksc
  meta_write <- meta_write[complete.cases(meta_write),]
  
  keeps <- c('feature', as.character(
    meta_write$specimenID
  ))
  
  write.table(counts_write[,keeps],
              file = paste0('MAYO_', tis,'_Braak_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('MAYO_', tis,'_Braak_counts.txt'),
    name = paste0('MAYO ', tis, ' Sageseqr Input Braak Counts'),
    parentId=activity$properties$id ),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0('MAYO_', tis,'_Braak_counts.txt'))
  
  #Thal
  meta_write <- metadata_sageqr[as.character(metadata_sageqr$tissue) == tis,]
  meta_write$thal <- metadata[meta_write$specimenID,]$thal
  meta_write <- meta_write[complete.cases(meta_write),]
  
  keeps <- c('feature', as.character(
    meta_write$specimenID
  ))
  
  write.table(counts_write[,keeps],
              file = paste0('MAYO_', tis,'_Thal_counts.txt'),
              row.names = F,
              col.names = T,
              quote = F,
              sep = '\t'
  )
  
  ENRICH_OBJ <- synapser::synStore( synapser::File(
    path=paste0('MAYO_', tis,'_Thal_counts.txt'),
    name = paste0('MAYO ', tis, ' Sageseqr Input Thal Counts'),
    parentId=activity$properties$id ),
    used = counts_used,
    activityName = activityName,
    executed = thisFile,
    activityDescription = activityDescription
  )
  synapser::synSetAnnotations(ENRICH_OBJ, annotations = all.annotations.expression)
  file.remove(paste0('MAYO_', tis,'_Thal_counts.txt'))
}
