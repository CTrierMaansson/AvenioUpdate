#' Add a new AVENIO run to patient overview
#' 
#' This function takes the main list of NGS data and adds the information from
#' a new run to list. If the patient ID/CPR is not present in the main list
#' the patient is added to the list. If a patient ID/CPR is present in the list 
#' the new NGS data is added under that patient and all the BAM files for that
#' patient is manually checked for mutations. 
#' @importFrom dplyr select `%>%` filter
#' @importFrom BiocBaseUtils isScalarCharacter
#' @importFrom readxl read_xlsx
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
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
    `%ni%` <- Negate(`%in%`)
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
    unlisted_before <- do.call(rbind,master_list) %>% 
        dplyr::select(sample_index,Analysis.ID,Sample.ID) %>% 
        unique()
    message("Reading run information")
    ID <- gsub(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-",
        "",Directory)
    run_ID_short <- substr(ID,1,8)
    AVENIO_runs <- readxl::read_xlsx(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    runs_ID <- AVENIO_runs %>% 
        dplyr::filter(nchar(Run_ID) != 24)
    runs_project <- AVENIO_runs %>% 
        filter(grepl("_",Project))
    if(nrow(runs_project)>0){
        failed_Pojects <- paste(unique(runs_project$Project),collapse = ", ")
        stop("Error in //Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx\n",
             "The following projects contain '_' which is not allowed:\n",
             failed_Pojects,"\n",
             "Please fix before proceeding.")
    }
    if(nrow(runs_ID)>0){
        failed_ID <- paste(runs_ID$Run_ID,collapse = ", ")
        stop("Error in //Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx\n",
             "The following run IDs do not contain the mandatory 24 characters: ",
             failed_ID)
    }
    AVENIO_runs_select <- AVENIO_runs %>% 
        dplyr::filter(grepl(run_ID_short,Run_ID)) %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date))
    if(any(is.na(AVENIO_runs_select$Project))){
        na_Project <- AVENIO_runs_select %>% 
            dplyr::filter(is.na(Project))
        na_Project_samples <- paste(na_Project$Sample_name,collapse = ", ")
        warning(paste0("The following samples are lacking Project information",
                       " and will be excluded from the analysis:\n",
                       na_Project_samples))
    }
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
        na_name_samples <- paste(na_name$Sample_name,collapse = ", ")
        warning(paste0("The following samples are lacking project name information",
                       " and will be excluded from the analysis:\n",
                       na_name_samples))
    }
    if(any(is.na(AVENIO_runs_select$Material))){
        na_material <- AVENIO_runs_select %>% 
            dplyr::filter(is.na(Sample_date))
        na_material_samples <- paste(na_material$Material,collapse = ", ")
        warning(paste0("The following samples are lacking source material ",
                       "information and will be excluded from the analysis:\n",
                       na_material_samples))
    }
    if(any(is.na(AVENIO_runs_select$date_check))){
        na_date <- AVENIO_runs_select %>% 
            dplyr::filter(is.na(date_check))
        na_date_samples <- paste(na_date$Sample_name,collapse = ", ")
        warning(paste0("The following samples have lacking date information,",
                       " or the date is not formatted correctly (YYYY-MM-DD)",
                       " and will be excluded from the analysis:\n",
                       na_date_samples))
    }
    AVENIO_runs_select <- AVENIO_runs_select %>% 
        dplyr::filter(!is.na(CPR)) %>% 
        dplyr::filter(!is.na(date_check)) %>% 
        dplyr::select(-date_check) %>% 
        dplyr::filter(!is.na(Name_in_project)) %>% 
        dplyr::filter(!is.na(Project))
    if(nrow(AVENIO_runs_select) == 0){
        stop("None of the samples have all the relevant information.\nSee the warnings above for details")
    }
    message("Merging run information and patient information")
    samples <- add_samples(Directory,AVENIO_runs_select)
    if(!any(samples$sample_index %ni% unlisted_before$sample_index)){
        stop("All samples are already part of the dataset. Terminating")
    }
    n_patients_before <- length(master_list)
    n_runs_before <- nrow(unlisted_before)
    message("Creating list of data.frames for new samples")
    df_list <- create_df_list(samples,AVENIO_runs_select)
    message("Merging existing data with the new run and verifies mutations in BAM files")
    reanalyzed <- reanalyze_samples(master_list,df_list)
    if(is.list(reanalyzed)){
        if(length(reanalyzed) >= length(master_list)){
            n_patients_after <- length(reanalyzed)
            unlisted_after <- do.call(rbind,reanalyzed) %>% 
                dplyr::select(Analysis.Name,Sample.ID) %>% 
                unique()
            n_runs_after <- nrow(unlisted_after)
            message(paste0("Before the dataset consisted of ",
                         n_patients_before, 
                         " individuals and ",
                         n_runs_before, 
                         " samples analyzed"))
            message(paste0("Now the dataset consists of ",
                         n_patients_after, 
                         " individuals and ",
                         n_runs_after, 
                         " samples analyzed"))
            samples_export_df <- AVENIO_runs_select %>% 
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
                                          sample_index)) %>% 
                dplyr::left_join(samples,by = "sample_index") %>% 
                dplyr::select(sample_index,Project) %>%
                unique()
            samples_count_project <- samples_export_df %>% 
                dplyr::count(Project)
            message("The following projects have been updated with this many samples")
            print(samples_count_project)
            message("And the following samples have been added to the dataset")
            print(samples_export_df$sample_index)
            message("Saving updated list of patients")
            saveRDS(reanalyzed,file = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
        }
    }
    return(reanalyzed)
}