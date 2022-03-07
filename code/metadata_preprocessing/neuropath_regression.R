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
  ref <- read.csv(synapser::synGet('syn25792425')$path,header=T) %>%
    tibble::column_to_rownames('specimenID')
  meta$individualID <- ref[meta$specimenID,]$individualID
  
  #Filter to match counts
  counts <- data.table::fread(
    synapser::synGet('syn25808335')$path, 
    sep = '\t',
    header = T) %>%
    as.data.frame()
  
  
  counts <- counts[, c('feature',meta$specimenID)]
  return(list(meta = meta,
              counts = counts
        ))
  
}

# RosMap

#Braak
braak <- process( 
  syns = list(
    ACC = 'syn26528995',
    PCC = 'syn26529011',
    DLPFC = 'syn26529008' 
    ),
  counts = 'syn25808335',
  inds = 'syn25792425'
)

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
  thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath='code/metadata_preprocessing/neuropath_regression.R')
  
  for( file in c('meta','counts')){
    # Write File
    if(file == 'meta') {
      write.csv(
        object[[file]], 
        file = paste0(file_names[file], '.',file_handle[file]), 
        row.names = F, 
        col.names = T, 
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
      executed = , 
      activityDescription = activ[[file]][['description']]
    )
    #Copy annotations from existing objects
    synapser::synSetAnnotations(
      ENRICH_OBJ, 
      annotations = synapser::synGetAnnotations(annots[[file]])
    )
  }
}

synapse_push(
  object = braak,
  parentID = list(meta = 'syn27554188', counts = 'syn27554176'),
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


