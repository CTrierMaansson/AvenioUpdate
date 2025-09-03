#' @noRd
#' @importFrom dplyr arrange `%>%` filter mutate select bind_rows mutate_all
#' @importFrom stringr str_split_i
#' @importFrom readxl read_xlsx
reanalyze_samples <- function(master_list, df_list,synology_path){
    `%ni%` <- Negate(`%in%`)
    existing <- names(df_list)[names(df_list) %in% names(master_list)]
    if(length(existing)<1){
        df_list <- lapply(df_list, function(df){
            df[] <- lapply(df, as.character)
            return(df)
        })
        test_df <- do.call(dplyr::bind_rows,df_list) %>% 
            dplyr::bind_rows(master_sele_df) %>% 
            dplyr::arrange(sample_index,Gene) %>% 
            unique()
    }else{
        master_sele <- master_list[names(master_list) %in% existing]
        master_sele <- lapply(master_sele, function(df){
            df[] <- lapply(df, as.character)
            return(df)
        })
        master_sele_df <- do.call(dplyr::bind_rows,master_sele)
        df_list <- lapply(df_list, function(df){
            df[] <- lapply(df, as.character)
            return(df)
        })
        test_df <- do.call(dplyr::bind_rows,df_list) %>% 
            dplyr::bind_rows(master_sele_df) %>% 
            dplyr::arrange(sample_index,Gene) %>% 
            unique()
    }
    rownames(test_df) <- NULL
    validated_samples <- add_variants_from_bam_files(test_df,synology_path)
    Avenio_runs <- readxl::read_xlsx(
        paste0(synology_path,"AVENIO_runs.xlsx"),
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    Avenio_sele <- create_sample_index(Avenio_runs) %>% 
        dplyr::filter(sample_index %in% validated_samples$sample_index)
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
            certain_BC <- sample_BC %>% 
                filter(!grepl("BAM",Flags) | Variant.Depth > 2)
            uncertain_BC <- sample_BC %>% 
                filter(grepl("BAM",Flags) & Variant.Depth < 3)
            sample_other <- sample_df %>% 
                dplyr::filter(!grepl("_BC", sample_index))
            if(nrow(sample_other) > 0){
                certain_BC <- certain_BC %>% 
                    dplyr::mutate(identifier = paste0(
                        Gene,
                        "_",
                        Coding.Change,
                        "_",
                        Genomic.Position))
                uncertain_BC <- uncertain_BC %>% 
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
                                      ifelse(identifier %in% certain_BC$identifier,
                                             ifelse(grepl("BC_mut",Flags),
                                                    Flags,ifelse(nchar(Flags)>0,
                                                                 paste0(Flags,", BC_mut"),
                                                                 paste0("BC_mut"))),
                                             Flags)) %>%
                    dplyr::mutate(Flags = 
                                      ifelse(identifier %in% uncertain_BC$identifier,
                                             ifelse(grepl("BC_mut",Flags),
                                                    Flags,ifelse(nchar(Flags)>0,
                                                                 paste0(Flags,", uncertain_BC_mut"),
                                                                 paste0("uncertain_BC_mut"))),
                                             Flags)) %>%
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