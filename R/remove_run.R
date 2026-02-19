#' Remove all results for a specific sequencing run 
#' 
#' Remove results for a run_ID in cases of e.g., wrong entries in 
#' 'Avenio_runs.xlsx', errors during sequencing or the results are misleading
#' @importFrom BiocBaseUtils isScalarCharacter isScalarLogical isCharacter
#' @importFrom dplyr `%>%` filter as_tibble
#' @param run `Character` or `Character` `vector` of run_IDs you want
#'  removed from the results
#' @param master_list The `list` of `data.frames` which contains all information on
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/"
#' @param save `Boolean` Indicating whether the output without the specified 
#'  `run` should be automatically saved as 
#'  "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds". 
#'  Default = `FALSE`
#' @param save_as Full path to the file you want to save the output as. 
#'  Default = paste0(synology_path,"AVENIO_results_patients.rds")
#' @return A `list` of `data.frames` with similar structure to `master_list`.
#'  but with the run_ID defined with `run` removed from the results
#'  
#' @examples
#' master <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' #Removing a single run from the results automatically
#' rm_single_run <- "EXAMPLE_RUN123"
#' remove_run(run = rm_single_run,
#'            master_list = master,
#'            save = TRUE)
#' 
#' #Removing several runs but without saving the results
#' rm_several_runs <- c("EXAMPLE_RUN123",
#'                      "EXAMPLE_RUN456",
#'                      "EXAMPLE_RUN789")
#' remove_run(run = rm_several_runs,
#'            master_list = master,
#'            save = FALSE)
#' @export
remove_run <- function(run, 
                       master_list,
                       synology_path = "//Synology_m1/Synology_folder/AVENIO/",
                       save = FALSE,
                       save_as = paste0(synology_path,
                                        "AVENIO_results_patients.rds")){
    `%ni%` <- Negate(`%in%`)
    if (!is.list(master_list)) {
        stop("master_list has to be a list")
    }
    for(i in 1:length(master_list)){
        if(!is.data.frame(master_list[[i]])){
            stop(paste0("Entry ",i," in the master_list is not a data.frame"))
        }
    }
    if (!isScalarCharacter(synology_path)) {
        stop("synology_path has to be a character")
    }
    if (!dir.exists(synology_path)) {
        stop("The path entered as synology_path does not exist")
    }
    nchar_path <- nchar(synology_path)
    if(substr(synology_path,nchar_path,nchar_path) != "/"){
        stop("The synology_path has to end with a '/'")
    }
    if (!isScalarLogical(save)) {
        stop("save has to be a TRUE or FALSE")
    }
    if (!isCharacter(run)) {
        stop("run has to be a character")
    }
    nchar_saveas <- nchar(save_as)
    if(substr(save_as,1,1) != "/"){
        stop("The save_as has to start with a '/'")
    }
    if(substr(save_as,nchar_saveas-3,nchar_saveas) != ".rds"){
        stop("The save_as has to end with .rds")
    }
    df <- unlist_frames(master_list = master_list)
    initial_runs <- unique(df$Analysis.ID)
    if(any(run %ni% initial_runs)){
        wrong_runs <- paste(unique(run[run %ni% initial_runs]),collapse = ", ")
        stop("Error in run\n",
             "The following runs that you want to remove\n",
             "are not present in master_list:\n",
             wrong_runs,"\n")
    }
    message(paste0("Removing\n",paste(run,collapse = ", "),
                   "\nfrom master_list\n\n"))
    df_filtered <- df %>% 
        filter(Analysis.ID %ni% run)
    df_removed <- df %>% 
        filter(Analysis.ID %in% run) %>% 
        select(CPR,sample_index, Analysis.ID, Flags, Gene,Coding.Change) %>% 
        unique() %>% 
        as_tibble()
    message("The following results have been removed")
    print(df_removed, n = nrow(df_removed))
    message("Recreating the list of data.frames with results")
    unique_patients <- unique(df_filtered$CPR) 
    df_list <- list()
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- df_filtered %>%
            dplyr::filter(CPR == CPR_patient) %>%
            select(-CPR)
        df_list[[i]] <- patient_df
        names(df_list)[i] <- CPR_patient
    }
    message("Checking correct list of results")
    df_check <- unlist_frames(master_list = df_list)
    index_check <- all.equal(df_filtered$sample_index,df_check$sample_index)
    if(!isScalarLogical(index_check)){
        stop(paste0("Error in recreated list of data.frames\n",
                    "This is the output of\n",
                    "all.equal(df_filtered$sample_index,df_check$sample_index)\n",
                    index_check[1], "\n",
                    index_check[2]))
    }
    else{
        if(index_check){
            message("Check successful")
        }
    }
    if(!save){
        n <- save_as_func(save_as_path = save_as)
        if(n == "1"){
            file_exist_func(save_as_path = save_as,
                            obj = df_list)
        }
        if(n == "2"){
            message(paste0("Terminating code without saving the output ",
                           "as an .rds file\n",
                           "But the filtered results are returned by this function"))
        }
        
    }
    if(save){
        file_exist_func(save_as_path = save_as,
                        obj = df_list)
    }
    return(df_list)
}