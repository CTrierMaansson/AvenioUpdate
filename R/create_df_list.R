#' @noRd
#' @importFrom readxl read_xlsx
#' @importFrom dplyr filter `%>%`
create_df_list <- function(sample_df){
    samples <- unique(sample_df$Sample.ID)
    AVENIO_runs <- readxl::read_xlsx(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    AVENIO_runs_select <- AVENIO_runs %>% 
        dplyr::filter(Sample_name %in% samples)
    unique_patients <- unique(AVENIO_runs_select$CPR)
    df_list <- list()
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- AVENIO_runs_select %>% 
            dplyr::filter(CPR == CPR_patient)
        patient_variants <- sample_df %>% 
            dplyr::filter(Sample.ID %in% patient_df$Sample_name)
        df_list[[i]] <- patient_variants 
        names(df_list)[i] <- CPR_patient
    }
    return(df_list)
}