#' Add a new AVENIO run to patient overview
#' 
#' This function takes the main list of NGS data and adds the information from
#' a new run to list. If the patient ID/CPR is not present in the main list
#' the patient is added to the list. If a patient ID/CPR is present in the list 
#' the new NGS data is added under that patient and all the BAM files for that
#' patient is manually checked for mutations. 
#' @importFrom dplyr select `%>%`
#' @importFrom BiocBaseUtils isScalarCharacter
#' @param master_list The `list` of `data.frames` which contains all information on
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @param  Directory `Character` string with the full path to the directory on
#'  the Synology where the data is located.
#' @return A `list` of `data.frames` with similar structure to `master_list`. 
#'  During the analysis the number of patients and NGS runs present in the 
#'  dataset before and after the NGS run is printed. If the list successfully 
#'  has been updated the `.rds` file on the Synology containing the list is 
#'  automatically updated. 
#'  
#' @examples
#' master <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' test_path <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ANei_Hsi__lGmaJ1k_Ji_d3O"
#' add_run_to_list(master_list = master,
#'                 Directory = test_path)
#' @export
add_run_to_list <- function(master_list, Directory){
    if (!is.list(master_list)) {
        stop("df_list has to be a list")
    }
    for(i in 1:length(master_list)){
        if(!is.data.frame(master_list[[i]])){
            stop(paste0("Entry ",i," in the master_list is not a data.frame"))
        }
    }
    if (!isScalarCharacter(Directory)) {
        stop("Directory has to be a character")
    }
    if (!grepl("//Synology",Directory)) {
        stop("'//Synology' has to be the first characters in the Directory string")
    }
    if (!file.exists(Directory)) {
        stop("The path entered as Directory does not exist")
    }
    samples <- add_samples(Directory)
    n_patients_before <- length(master_list)
    unlisted_before <- do.call(rbind,master_list) %>% 
        dplyr::select(Analysis.Name,Sample.ID) %>% 
        unique()
    n_runs_before <- nrow(unlisted_before)
    df_list <- create_df_list(samples)
    reanalyzed <- reanalyze_samples(master_list,df_list)
    if(is.list(reanalyzed)){
        if(length(reanalyzed) >= length(master_list)){
            n_patients_after <- length(reanalyzed)
            unlisted_after <- do.call(rbind,reanalyzed) %>% 
                dplyr::select(Analysis.Name,Sample.ID) %>% 
                unique()
            n_runs_after <- nrow(unlisted_after)
            print("Before the dataset consisted of ",
                  n_patients_before, 
                  " patients and ",
                  n_runs_before, 
                  " samples analyzed")
            print("Now the dataset consists of ",
                  n_patients_after, 
                  " patients and ",
                  n_runs_after, 
                  " samples analyzed")
            saveRDS(reanalyzed,file = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
        }
    }
    return(reanalyzed)
}