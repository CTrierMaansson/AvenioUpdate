#' @noRd
#' @importFrom dplyr arrange `%>%` filter mutate select
#' @importFrom stringr str_split_i
#' @importFrom readxl read_xlsx
reanalyze_samples <- function(master_list, df_list){
    `%ni%` <- Negate(`%in%`)
    existing <- names(df_list)[names(df_list) %in% names(master_list)]
    if(length(existing)<1){
        test_df <- do.call(rbind,df_list) %>% 
            dplyr::arrange(sample_index,Gene) %>% 
            unique()
    }
    else{
        master_sele <- master_list[names(master_list) %in% existing]
        test_df <- do.call(rbind,append(master_sele,df_list)) %>% 
            dplyr::arrange(sample_index,Gene) %>% 
            unique()
    }
    
    rownames(test_df) <- NULL
    validated_samples <- add_variants_from_bam_files(test_df)
    Avenio_runs <- readxl::read_xlsx(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    Avenio_sele <- Avenio_runs %>% 
        dplyr::filter(Run_ID %in% validated_samples$Analysis.ID) %>% 
        dplyr::filter(Sample_name %in% validated_samples$Sample.ID)
    validated_list <- create_df_list(validated_samples,Avenio_sele)
    for(i in 1:length(validated_list)){
        sample_df <- validated_list[[i]]
        if(any(grepl("_BC",sample_df$sample_index))){
            patient_sample <- paste0(
                stringr::str_split_i(unique(
                    sample_df$sample_index)[1],"_",1),
                "_",
                stringr::str_split_i(unique(
                    sample_df$sample_index)[1],"_",2))
            message(paste0("There is a BC sample for ",
                         patient_sample,
                         " and variants are marked accordingly"))
            sample_BC <- sample_df %>% 
                dplyr::filter(grepl("_BC", sample_index))
            sample_other <- sample_df %>% 
                dplyr::filter(!grepl("_BC", sample_index))
            if(nrow(sample_other) > 0){
                sample_BC <- sample_BC %>% 
                    dplyr::mutate(identifier = paste0(
                        Gene,
                        "_",
                        Coding.Change,
                        "_",
                        Genomic.Position))
                sample_other <- sample_other %>% 
                    dplyr::mutate(identifier = paste0(
                        Gene,
                        "_",
                        Coding.Change,
                        "_",
                        Genomic.Position)) %>% 
                    dplyr::mutate(Flags = ifelse(is.na(Flags),"",Flags)) %>% 
                    dplyr::mutate(Flags = 
                                      ifelse(identifier %in% sample_BC$identifier,
                                             ifelse(grepl("BC_mut",Flags),
                                                    Flags,ifelse(nchar(Flags)>0,
                                                                 paste0(Flags,", BC_mut"),
                                                                 paste0("BC_mut"))),
                                             Flags)) %>%  
                    dplyr::select(-identifier)
                sample_BC <- sample_BC %>% 
                    dplyr::select(-identifier)
                sample_df <- rbind(sample_other,sample_BC) %>% 
                    dplyr::arrange(sample_index,Gene)
            }
        }
        validated_list[[i]] <- sample_df
    }
    if(length(existing)>0){
        master_list <- master_list[names(master_list) %ni% names(validated_list)]
        master_list <- append(master_list,validated_list)
    }
    else{
        master_list <- append(master_list,validated_list)
    }
    master_list <- master_list[sort(names(master_list))]
    return(master_list)
}