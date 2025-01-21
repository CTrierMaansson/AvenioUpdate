#' @noRd
#' @importFrom utils unzip
#' @importFrom readxl read_xlsx
#' @importFrom dplyr filter select full_join left_join mutate relocate arrange `%>%`
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
add_samples <- function(Directory){
    `%ni%` <- Negate(`%in%`)
    Run_ID <- gsub(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-",
        "",Directory)
    zip_file <- paste0(Directory,"/","all_CSVs-",Run_ID,".zip")
    sample_file <- paste0(Directory,"/SampleMetrics.csv")
    run_ID_short <- substr(Run_ID,1,8)
    unzip(zipfile = zip_file,
          exdir = Directory)
    Filter_variants <- paste0(Directory,"/","VariantMetrics-UserFilterSet.csv")
    variants <- read.csv(Filter_variants)
    AVENIO_runs <- readxl::read_xlsx(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    AVENIO_runs_select <- AVENIO_runs %>% 
        dplyr::filter(grepl(run_ID_short,Run_ID)) %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date))
    if(any(is.na(AVENIO_runs_select$date_check))){
        na_date <- AVENIO_runs_select %>% 
            dplyr::filter(is.na(date_check))
        na_date_samples <- paste(na_date$Sample_name,collapse = ", ")
        warning(paste0("The following samples have lacking date information,",
                       " or the date is not formatted correctly (YYYY-MM-DD)",
                       " and will be excluded from the analysis:\n",
                       na_date_samples))
    }
    variants_select <- variants %>% 
        dplyr::filter(Sample.ID %in% AVENIO_runs_select$Sample_name)
    samples_with_variants <- unique(variants_select$Sample.ID)
    samples_of_interest <- unique(AVENIO_runs_select$Sample_name)
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
        if(any(is.na(AVENIO_runs_select$CPR))){
            na_CPR <- AVENIO_runs_select %>% 
                dplyr::filter(is.na(CPR))
            na_CPR_samples <- paste(na_CPR$Sample_name,collapse = ", ")
            warning(paste0("The following samples are lacking CPR information",
                           " and will be excluded from the analysis:\n",
                           na_CPR_samples))
        }
        if(any(is.na(AVENIO_runs_select$Name_in_project))){
            na_name <- AVENIO_runs_select %>% 
                dplyr::filter(is.na(Name_in_project))
            na_name_samples <- paste(na_date$Sample_name,collapse = ", ")
            warning(paste0("The following samples are lacking project name information",
                           " and will be excluded from the analysis:\n",
                           na_name_samples))
        }
        if(any(is.na(AVENIO_runs_select$Material))){
            na_material <- AVENIO_runs_select %>% 
                dplyr::filter(is.na(Sample_date))
            na_material_samples <- paste(na_date$Material,collapse = ", ")
            warning(paste0("The following samples are lacking source material ",
                           "information and will be excluded from the analysis:\n",
                           na_material_samples))
        }
        AVENIO_runs_select <- AVENIO_runs_select %>%
            dplyr::filter(!is.na(CPR)) %>% 
            dplyr::filter(!is.na(date_check)) %>% 
            dplyr::select(-date_check) %>% 
            dplyr::filter(!is.na(Name_in_project)) %>% 
            dplyr::mutate(sample_index = paste0(Project,"_",
                                         Name_in_project,"_",
                                         substr(stringr::str_split_i(
                                             as.character(Sample_date),"-",1),
                                             3,4),
                                         stringr::str_split_i(
                                             as.character(Sample_date),"-",2),
                                         stringr::str_split_i(
                                             as.character(Sample_date),"-",3)
                                         )
                   ) %>% 
            dplyr::mutate(
                sample_index = ifelse(Material != "cfDNA",
                                      paste0(sample_index,"_",Material),
                                      sample_index))
    }
    AVENIO_runs_select_merge <- AVENIO_runs_select %>% 
        dplyr::select(Sample_name,sample_index)
    combined_samples <- combined_samples %>% 
        dplyr::filter(Sample.ID %in% AVENIO_runs_select$Sample_name) %>% 
        dplyr::left_join(AVENIO_runs_select_merge,
                         by = c("Sample.ID"="Sample_name")) %>% 
        dplyr::relocate(sample_index,.before = Sample.ID) %>% 
        dplyr::arrange(sample_index,Gene)
    return(combined_samples)
}