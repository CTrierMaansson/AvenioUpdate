#' Explores the .xlsx file with the run information
#' 
#' This function explores the AVENIO_runs.xlsx file and performs various tests
#' on the information and links the information to the list of mutations
#' for each patient 
#' @importFrom dplyr filter `%>%` mutate select group_by ungroup count arrange
#' @importFrom BiocBaseUtils isScalarCharacter isScalarLogical
#' @importFrom readxl read_xlsx
#' @importFrom lubridate ymd
#' @param  Info `Character` string with the name of the information that should.
#' be extracted. Possible options can be viewed with explore_AVENIO_runs_Info().
#' If `NULL` (default) all information is returned as a `list`.
#' @param silent `Boolean` determining information messages should be displayed 
#' If `FALSE` (default) messages are displayed. If `TRUE` the messages are
#' not printed. 
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/"  
#' @return A `list` with different information on the AVENIO_runs.xslx 
#' (`Info` = `NULL`) or the specific object as determined by `Info`.
#'  
#' @examples
#' results <- explore_AVENIO_runs()
#' results$Total_entries
#' explore_AVENIO_runs(Info = "Required")
#' @export
explore_AVENIO_runs <- function(Info = NULL, silent = FALSE,
                                synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    `%ni%` <- Negate(`%in%`)
    if(!is.null(Info)){
        if (!isScalarCharacter(Info)) {
            stop("Info has to be a character or NULL")
        }
    }
    if (!isScalarLogical(silent)) {
        stop("silent has to be a TRUE or FALSE")
    }
    if(!silent){
        message("Reading AVENIO_runs.xlsx")
    }
    df <- readxl::read_xlsx(
        paste0(synology_path,"AVENIO_runs.xlsx"),
        col_types = c(rep("guess",4),"date",rep("guess",6))) %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date))
    if(!silent){
        message("Exploring complete and incomplete information")
    }
    total_entries <- nrow(df)
    runs_ID_incomplete <- df %>% 
        dplyr::filter(nchar(Run_ID) != 24)
    if(nrow(runs_ID_incomplete) == 0){
        runs_ID_incomplete <- NULL
    }
    required_information <- df %>% 
        dplyr::filter(!is.na(CPR)) %>% 
        dplyr::filter(!is.na(Name_in_project)) %>% 
        dplyr::filter(!is.na(Project)) %>% 
        dplyr::filter(!is.na(date_check)) %>% 
        dplyr::filter(!is.na(Material)) %>% 
        dplyr::filter(!is.na(Run_ID)) %>% 
        dplyr::filter(!is.na(Run_name)) %>% 
        dplyr::filter(!is.na(Sample_name)) %>% 
        dplyr::select(CPR,Name_in_project,
                      Project,Sample_date,
                      Run_name,Run_ID,
                      Sample_name,Material,
                      Sample_note)
    incomplete_dates <- df %>% 
        dplyr::filter(is.na(date_check))
    if(nrow(incomplete_dates) == 0){
        incomplete_dates <- NULL
    }
    incomplete_names <- df %>% 
        dplyr::filter(is.na(Name_in_project) | is.na(Project) | is.na(Sample_name))
    if(nrow(incomplete_names) == 0){
        incomplete_names <- NULL
    }
    incomplete_material <- df %>% 
        dplyr::filter(is.na(Material))
    if(nrow(incomplete_material) == 0){
        incomplete_material <- NULL
    }
    if(!silent){
        message("Extracting stats")
    }
    mat_stats <- required_information %>% 
        dplyr::select(Material,CPR) %>%
        dplyr::group_by(Material) %>% 
        dplyr::count() %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(desc(n))
    note_stats <- required_information %>% 
        dplyr::select(Sample_note,CPR) %>%
        dplyr::group_by(Sample_note) %>% 
        dplyr::count() %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(desc(n))
    if(!silent){
        message("Reading AVENIO_results_patients.rds")
    }
    master_list <- readRDS(paste0(synology_path,"AVENIO_results_patients.rds"))
    if(!silent){
        message("Investigating AVENIO_results_patients.rds")
    }
    analyses <- AvenioUpdate::included_analyses(master_list)
    non_included_IDs <- df %>%
        dplyr::filter(!is.na(Run_ID)) %>% 
        dplyr::filter(Run_ID %ni% analyses[[1]]$Analysis.ID) %>% 
        dplyr::select(Project,Run_name, Run_ID) %>% 
        unique()
    if(nrow(non_included_IDs) == 0){
        non_included_IDs <- NULL
    }
    non_CPRs <- df %>%
        dplyr::filter(!is.na(CPR)) %>% 
        dplyr::filter(CPR %ni% names(master_list)) %>% 
        dplyr::select(CPR,Project) %>% 
        unique()
    if(nrow(non_CPRs) == 0){
        non_CPRs <- NULL
    }
    if(!silent){
        message("Formatting output")
    }
    result = list(Total_entries = total_entries,
                  Complete_entries = nrow(required_information),
                  Required = required_information,
                  Material_stats = mat_stats,
                  Time_stats = note_stats,
                  Unincluded_analyses = non_included_IDs,
                  Unincluded_CPRs = non_CPRs,
                  Incomplete_IDs = runs_ID_incomplete,
                  Incomplete_dates = incomplete_dates,
                  Incomplete_names = incomplete_names,
                  Incomplete_material = incomplete_material)
    
    if(!is.null(Info)){
        if(Info %ni% names(result)){
            stop("Info has to a a viable name of the output list. Relevant entries can be viewed with explore_AVENIO_runs_Info()")
        }
        if(!silent){
            message(paste0("Outputting: ",Info))
        }
        output <- result[[Info]]
    }
    else{
        output <- result
    }
    return(output)
}