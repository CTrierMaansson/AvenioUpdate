#' @noRd
#' @importFrom dplyr filter `%>%`
create_df_list <- function(sample_df,runs){
    samples <- unique(sample_df$Sample.ID)
    unique_patients <- unique(runs$CPR)
    df_list <- list()
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- runs %>% 
            dplyr::filter(CPR == CPR_patient)
        patient_variants <- sample_df %>% 
            dplyr::filter(Sample.ID %in% patient_df$Sample_name)
        df_list[[i]] <- patient_variants 
        names(df_list)[i] <- CPR_patient
    }
    return(df_list)
}