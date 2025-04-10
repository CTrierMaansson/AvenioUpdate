#' @noRd
#' @importFrom dplyr filter `%>%` mutate select
#' @importFrom stringr str_split_i 
create_df_list <- function(sample_df,runs){
    samples <- unique(sample_df$sample_index)
    unique_patients <- unique(runs$CPR)
    df_list <- list()
    runs <- runs %>% 
        dplyr::mutate(sample_index = paste0(Project,"_",
                                            Name_in_project,"_",
                                            substr(stringr::str_split_i(
                                                as.character(Sample_date),"-",1),
                                                3,4),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",2),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",3))) %>% 
        dplyr::mutate(
            sample_index = ifelse(Material != "cfDNA",
                                  paste0(sample_index,"_",Material),
                                  sample_index))
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- runs %>% 
            dplyr::filter(CPR == CPR_patient)
        patient_variants <- sample_df %>% 
            dplyr::filter(sample_index %in% patient_df$sample_index) %>% 
            dplyr::select(-sample_index)
        df_list[[i]] <- patient_variants 
        names(df_list)[i] <- CPR_patient
    }
    return(df_list)
}