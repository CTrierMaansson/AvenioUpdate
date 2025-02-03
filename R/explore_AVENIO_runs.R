#' Explores the .xlsx file with the run information
#' 
#' This function explores the AVENIO_runs.xlsx file and performs various tests
#' on the information and links the information to the list of mutations
#' for each patient 
#' @importFrom dplyr filter `%>%` mutate
#' @importFrom BiocBaseUtils isScalarCharacter
#' @importFrom readxl read_xlsx
#' @param  Info `Character` string with the
#' @return A `list` with different information on the AVENIO_runs.xslx
#'  
#' @examples
#' results <- explore_AVENIO_runs()
#' results[["Total_entries"]]
#' explore_AVENIO_runs(Info = "Required")
#' @export
explore_AVENIO_runs <- function(Info = NULL){
    `%ni%` <- Negate(`%in%`)
    if(!is.null(Info)){
        if (!isScalarCharacter(Info)) {
            stop("Info has to be a character or NULL")
        }
    }
    df <- readxl::read_xlsx(
        "//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
        col_types = c(rep("guess",4),"date",rep("guess",6))) %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date))
    total_entries <- nrow(df)
    runs_ID_incomplete <- df %>% 
        dplyr::filter(nchar(Run_ID) != 24)
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
                      Sample_name,Material)
    
    result = list(Total_entries = total_entries,
                  Required_entries = nrow(required_information),
                  Required = required_information,
                  Incomplete_IDs = runs_ID_incomplete)
    
    if(!is.null(Info)){
        if(Info %ni% names(result)){
            stop("Info has to a a viable name of the output list. Relevant entries can be viewed with XYZ")
        }
        output <- result[[Info]]
    }
    else{
        output <- result
    }
    return(output)
}
