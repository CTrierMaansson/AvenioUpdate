#' @noRd
#' @importFrom utils unzip
#' @importFrom readxl read_xlsx
#' @importFrom dplyr filter select full_join left_join mutate relocate arrange `%>%`
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
add_samples <- function(Directory,runs,synology_path){
    `%ni%` <- Negate(`%in%`)
    bname <- basename(Directory)
    Analysis_type <- unlist(lapply(strsplit(bname,"-"), "[[", 1))
    if(Analysis_type %ni% c("Tissue","Plasma")){
        stop(paste0("The final directory entered as Directory has to ",
                    "start with 'Plasma-' or Tissue-"))
    }
    Run_ID <- gsub(
        paste0(synology_path,"AVENIO_results/Plasma-"),
        "",Directory)
    if(Analysis_type == "Plasma"){
        Run_ID <- gsub(
            paste0(synology_path,"AVENIO_results/Plasma-"),
            "",Directory)
    }    
    if(Analysis_type == "Tissue"){
        Run_ID <- gsub(
            paste0(synology_path,"AVENIO_results/Tissue-"),
            "",Directory)
    }
    zip_file <- paste0(Directory,"/","all_CSVs-",Run_ID,".zip")
    sample_file <- paste0(Directory,"/SampleMetrics.csv")
    run_ID_short <- substr(Run_ID,1,8)
    unzip(zipfile = zip_file,
          exdir = Directory)
    Filter_variants <- paste0(Directory,"/","VariantMetrics-UserFilterSet.csv")
    variants <- read.csv(Filter_variants)
    AVENIO_runs_select <- runs
    variants_select <- variants %>% 
        dplyr::filter(Sample.ID %in% AVENIO_runs_select$Sample_name)
    samples_with_variants <- unique(variants_select$Sample.ID)
    samples_of_interest <- unique(AVENIO_runs_select$Sample_name)
    sample_metrics <- read.csv(sample_file)
    if(any(samples_of_interest %ni% unique(sample_metrics$Sample.ID))){
        error_sample <- samples_of_interest[samples_of_interest %ni% 
                                                unique(sample_metrics$Sample.ID)]
        error_samples <- paste(error_sample,collapse = ", ")
        stop(paste0("Sample(s): ",
                    error_samples,
                    " have an incorrect name and the sample ID cannot be found"))
    }
    if(length(samples_with_variants) != length(samples_of_interest)){
        sample_metrics <- read.csv(sample_file) %>% 
            filter(Sample.ID %in% samples_of_interest) %>% 
            filter(Sample.ID %ni% samples_with_variants)
        sample_metrics <- sample_metrics %>% 
            dplyr::select(
                colnames(sample_metrics)[colnames(sample_metrics) %in% 
                                             colnames(variants_select)])
        sele_cols <- colnames(variants_select)[colnames(variants_select) %in%
                                                   colnames(sample_metrics)]
        for(i in 1:length(sele_cols)){
            if(class(sample_metrics[[sele_cols[i]]]) != 
               class(variants_select[[sele_cols[i]]])){
                classes <- c(class(sample_metrics[[sele_cols[i]]]),
                             class(variants_select[[sele_cols[i]]]))
                if("character" %in% classes){
                    sample_metrics[[sele_cols[i]]] <- as.character(
                        sample_metrics[[sele_cols[i]]])
                    variants_select[[sele_cols[i]]] <- as.character(
                        variants_select[[sele_cols[i]]])
                }
                if("numeric" %in% classes){
                    sample_metrics[[sele_cols[i]]] <- as.numeric(
                        sample_metrics[[sele_cols[i]]])
                    variants_select[[sele_cols[i]]] <- as.numeric(
                        variants_select[[sele_cols[i]]])
                }
            }
        }
        combined_samples <- dplyr::full_join(variants_select,
                                             sample_metrics,
                                             copy = TRUE)
    }
    else{
        combined_samples <- variants_select
    }
    if(!"sample_index" %in% colnames(AVENIO_runs_select)){
        AVENIO_runs_select <- create_sample_index(AVENIO_runs_select) 
    }
    AVENIO_runs_select_merge <- AVENIO_runs_select %>% 
        dplyr::select(Sample_name,sample_index) %>% 
        dplyr::mutate(Sample_name = as.character(Sample_name))
    combined_samples <- combined_samples %>% 
        dplyr::mutate(Sample.ID = as.character(Sample.ID)) %>% 
        dplyr::filter(Sample.ID %in% AVENIO_runs_select$Sample_name) %>%
        dplyr::left_join(AVENIO_runs_select_merge,
                         by = c("Sample.ID"="Sample_name")) %>% 
        dplyr::relocate(sample_index,.before = Sample.ID) %>% 
        dplyr::arrange(sample_index,Gene)
    return(combined_samples)
}