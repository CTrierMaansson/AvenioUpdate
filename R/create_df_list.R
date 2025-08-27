#' @noRd
#' @importFrom dplyr filter `%>%` mutate
#' @importFrom stringr str_split_i 
create_df_list <- function(sample_df,runs){
    samples <- unique(sample_df$sample_index)
    unique_patients <- unique(runs$CPR)
    df_list <- list()
    runs <- create_sample_index(runs)
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- runs %>% 
            dplyr::filter(CPR == CPR_patient)
        patient_variants <- sample_df %>% 
            dplyr::filter(sample_index %in% patient_df$sample_index)
        df_list[[i]] <- patient_variants 
        names(df_list)[i] <- CPR_patient
    }
    return(df_list)
}