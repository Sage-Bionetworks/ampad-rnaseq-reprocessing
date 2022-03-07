synapser::synLogin()
library(dplyr)

#' @param syns a named list of tissues and synIDs
#' @param counts Syn ID of Counts
#' @param inds Metadata with Indvidual IDs
process <- function( syns, counts, inds ){
  # Import Meta data and merge with tissue
  imported <- list()
  for(syn in names(syns)) { 
    imported[[syn]] <- read.csv(synapser::synGet(as.character(syns[syn]))$path, header = T)  
    imported[[syn]]$tissue <- syn
  }
  meta <- do.call( rbind, imported)
  
  #Annotate with indv. ID
  ref <- read.csv(synapser::synGet(inds)$path,header=T) %>%
    tibble::column_to_rownames('specimenID')
  meta$individualID <- ref[meta$specimenID,]$individualID
  
  #Filter to match counts
  counts <- data.table::fread(
    synapser::synGet(counts)$path, 
    sep = '\t',
    header = T) %>%
    as.data.frame()
  
  counts <- counts[, c('feature',meta$specimenID)]
  return(list(meta = meta,
              counts = counts
        ))
  
}

#' @param object a named list of `meta` and `counts`
#' @param parentID a `meta` and `counts` named list of folder parent IDs to push synapse meta and counts
#' @param used_syns the names of syn IDs to use for provenance
#' @param annots a `meta` and `counts` named list of synobjects to scrub annotations
#' @param file_names a `meta` and `counts` named list of names to
#' @param git_info a named list of `repo` and `file` that has the git repo and file name values
#' @param git_info a named list of `repo` and `file` that has the 

object <- braak
parentID <- list(meta = 'syn27554188',
                 counts = 'syn27554176')
used_syns <- c('syn26528995', 'syn26529011', 'syn26529008', 'syn25808335', 'syn25792425')
file_names <- list(meta = 'RosMap_Braak_Metadata',
                   counts = 'RosMap_Braak_Counts')
file_handle <- list(meta = 'csv',
                   counts = 'tsv')

annots <- list(meta = 'syn25808375',
               counts = 'syn25808335')
git_info <- list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                  file = 'code/metadata_preprocessing/neuropath_regression.R',
                  branch = 'main')
activ <- list(
    meta = list( name ='RosMap Braak Metadata',
                 description = 'RosMap metadata samples with Braak annotated'),
    counts = list( name ='RosMap Braak Counts',
                 description = 'RosMap Count Data matched to samples with Braak annotated'))
synapse_push <- function( object, parentID, used_syns, file_handle, file_names, annots, git_info, activ){
  
  #Github Info Pull
  thisRepo <- githubr::getRepo(repository = git_info[['repo']], ref="branch", refName=git_info[['branch']])
  thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath=git_info[['file']])
  
  for( file in c('meta','counts')){
    # Write File
    if(file == 'meta') {
      write.csv(
        object[[file]], 
        file = paste0(file_names[file], '.',file_handle[file]), 
        row.names = F, 
        quote = F)
    }else{
      write.table(
        object[[file]], 
        file = paste0(file_names[file], '.',file_handle[file]), 
        sep ='\t',  
        row.names = F, 
        col.names = T, 
        quote = F)
    }
    
    #Push File
    ENRICH_OBJ <- synapser::synStore( synapser::File( 
      path = paste0(file_names[file], '.',file_handle[file]), 
      name = file_names[[file]],
      parentId = parentID[[file]]),
      used = used_syns, 
      activityName = activ[[file]][['name']], 
      executed = thisFile, 
      activityDescription = activ[[file]][['description']]
    )
    #Copy annotations from existing objects
    synapser::synSetAnnotations(
      ENRICH_OBJ, 
      annotations = synapser::synGetAnnotations(annots[[file]])
    )
  }
}


# RosMap

## Braak
  braak <- process( 
    syns = list(
      ACC = 'syn26528995',
      PCC = 'syn26529011',
      DLPFC = 'syn26529008' 
    ),
    counts = 'syn25808335',
    inds = 'syn25792425'
  )
  
  synapse_push(
    object = braak,
    parentID = list(meta = 'syn27554311', counts = 'syn27554261'),
    used_syns = c('syn26528995', 'syn26529011', 'syn26529008', 'syn25808335', 'syn25792425'),
    file_names <- list(meta = 'RosMap_Braak_Metadata', counts = 'RosMap_Braak_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn25808375', counts = 'syn25808335'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                      file = 'code/metadata_preprocessing/neuropath_regression.R',
                      branch = 'main'),
    activ = list(
      meta = list( name ='RosMap Braak Metadata',
                   description = 'RosMap metadata samples with Braak annotated'),
      counts = list( name ='RosMap Braak Counts',
                     description = 'RosMap Count Data matched to samples with Braak annotated'))
  )

## CERAD
  cerad <- process( 
    syns = list(
      ACC = 'syn26529000',
      PCC = 'syn26529012',
      DLPFC = 'syn26529009' 
    ),
    counts = 'syn25808335',
    inds = 'syn25792425'
  )
  
  synapse_push(
    object = cerad,
    parentID = list(meta = 'syn27554311', counts = 'syn27554261'),
    used_syns = c('syn26529000', 'syn26529012', 'syn26529009', 'syn25808335', 'syn25792425'),
    file_names <- list(meta = 'RosMap_Cerad_Metadata', counts = 'RosMap_Cerad_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn25808375', counts = 'syn25808335'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='RosMap Cerad Metadata',
                   description = 'RosMap metadata samples with Cerad annotated'),
      counts = list( name ='RosMap Cerad Counts',
                     description = 'RosMap Count Data matched to samples with Cerad annotated'))
  )

## CogDX
  cogdx <- process( 
    syns = list(
      ACC = 'syn26529004',
      PCC = 'syn26529013',
      DLPFC = 'syn26529010' 
    ),
    counts = 'syn25808335',
    inds = 'syn25792425'
  )
  
  synapse_push(
    object = cogdx,
    parentID = list(meta = 'syn27554311', counts = 'syn27554261'),
    used_syns = c('syn26529004', 'syn26529013', 'syn26529010', 'syn25808335', 'syn25792425'),
    file_names <- list(meta = 'RosMap_CogDX_Metadata', counts = 'RosMap_CogDX_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn25808375', counts = 'syn25808335'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='RosMap CogDX Metadata',
                   description = 'RosMap metadata samples with CogDX annotated'),
      counts = list( name ='RosMap CogDX Counts',
                     description = 'RosMap Count Data matched to samples with CogDX annotated'))
  )

# Mayo
  braak <- process( 
    syns = list(
      TCX = 'syn26529080',
      CBE = 'syn26529086'
    ),
    counts = 'syn25835094',
    inds = 'syn26958125'
  )
  
  synapse_push(
    object = braak,
    parentID = list(meta = 'syn27554322', counts = 'syn27554273'),
    used_syns = c('syn26529080', 'syn26529086', 'syn25808335', 'syn25792425'),
    file_names <- list(meta = 'Mayo_Braak_Metadata', counts = 'Mayo_Braak_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn25835093', counts = 'syn25835094'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='Mayo Braak Metadata',
                   description = 'Mayo metadata samples with Braak annotated'),
      counts = list( name ='Mayo Braak Counts',
                     description = 'Mayo Count Data matched to samples with Braak annotated'))
  )
# Thal  
  thal <- process( 
    syns = list(
      TCX = 'syn26529081',
      CBE = 'syn26529087'
    ),
    counts = 'syn25835094',
    inds = 'syn26958125'
  )
  
  synapse_push(
    object = thal,
    parentID = list(meta = 'syn27554322', counts = 'syn27554273'),
    used_syns = c('syn26529081', 'syn26529087', 'syn25808335', 'syn25792425'),
    file_names <- list(meta = 'Mayo_Thal_Metadata', counts = 'Mayo_Thal_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn25835093', counts = 'syn25835094'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
      activ = list(
      meta = list( name ='Mayo Thal Metadata',
                   description = 'Mayo metadata samples with Braak annotated'),
      counts = list( name ='Mayo Thal Counts',
                     description = 'Mayo Count Data matched to samples with Thal annotated'))
  )
  
# MSBB
  # Braak  
  braak <- process( 
    syns = list(
      PHG = 'syn26525353',
      STG = 'syn26525356',
      IFG = 'syn26525347',
      FP = 'syn26525337'
    ),
    counts = 'syn25872241',
    inds = 'syn26525356'
  )
  
  synapse_push(
    object = braak,
    parentID = list(meta = 'syn27554333', counts = 'syn27554296'),
    used_syns = c('syn26525353', 'syn26525356', 'syn26525347', 'syn26525337', 'syn25872241', 'syn26525356'),
    file_names <- list(meta = 'MSBB_Braak_Metadata', counts = 'MSBB_Braak_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn26436062', counts = 'syn25872241'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='MSBB Braak Metadata',
                   description = 'MSBB metadata samples with Braak annotated'),
      counts = list( name ='MSBB Braak Counts',
                     description = 'MSBB Count Data matched to samples with Braak annotated'))
  )
  
  # Cerad  
  cerad <- process( 
    syns = list(
      PHG = 'syn26525354',
      STG = 'syn26525357',
      IFG = 'syn26525348',
      FP = 'syn26525341'
    ),
    counts = 'syn25872241',
    inds = 'syn26525356'
  )
  
  synapse_push(
    object = cerad,
    parentID = list(meta = 'syn27554333', counts = 'syn27554296'),
    used_syns = c('', '', '', '', 'syn25872241', 'syn26525356'),
    file_names <- list(meta = 'MSBB_Cerad_Metadata', counts = 'MSBB_Cerad_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn26436062', counts = 'syn25872241'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='MSBB Cerad Metadata',
                   description = 'MSBB metadata samples with Cerad annotated'),
      counts = list( name ='MSBB Cerad Counts',
                     description = 'MSBB Count Data matched to samples with Cerad annotated'))
  ) 
  
  # CogDx  
  cogdx <- process( 
    syns = list(
      PHG = 'syn26525355',
      STG = 'syn26525358',
      IFG = 'syn26525352',
      FP = 'syn26525346'
    ),
    counts = 'syn25872241',
    inds = 'syn26525356'
  )
  
  synapse_push(
    object = cogdx,
    parentID = list(meta = 'syn27554333', counts = 'syn27554296'),
    used_syns = c('syn26525355', 'syn26525358', 'syn26525352', 'syn26525346', 'syn25872241', 'syn26525356'),
    file_names <- list(meta = 'MSBB_CogDx_Metadata', counts = 'MSBB_CogDx_Counts'),
    file_handle = list(meta = 'csv',  counts = 'tsv'),
    annots = list(meta = 'syn26436062', counts = 'syn25872241'),
    git_info = list( repo ='Sage-Bionetworks/ampad-rnaseq-reprocessing',
                     file = 'code/metadata_preprocessing/neuropath_regression.R',
                     branch = 'main'),
    activ = list(
      meta = list( name ='MSBB CogDx Metadata',
                   description = 'MSBB metadata samples with CogDx annotated'),
      counts = list( name ='MSBB CogDx Counts',
                     description = 'MSBB Count Data matched to samples with CogDx annotated'))
  ) 
  
