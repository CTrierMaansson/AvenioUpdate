#' Highlights the parameters to be investigated from explore_AVENIO_runs()
#' 
#' This functions creates an overview of the different contents in the list
#' exported by explore_AVENIO_runs()
#' @importFrom dplyr `%>%` left_join
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/"
#' @return A `data.frame` with two variables containing the name and the
#' description of the entries in the explore_AVENIO_runs() output
#'  
#' @examples
#' explore_AVENIO_runs_Info()
#' @export
explore_AVENIO_runs_Info <- function(synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    if (!dir.exists(synology_path)) {
        stop("The path entered as synology_path does not exist")
    }
    nchar_path <- nchar(synology_path)
    if(substr(synology_path,nchar_path,nchar_path) != "/"){
        stop("The synology_path has to end with a '/'")
    }
    current_format <- explore_AVENIO_runs(silent = TRUE,
                                          synology_path = synology_path)
    df_current <- data.frame(Name = names(current_format))
    explanation <- data.frame(Name = c("Total_entries",
                                       "Complete_entries",
                                       "Required",
                                       "Material_stats",
                                       "Time_stats",
                                       "Unincluded_analyses",
                                       "Unincluded_CPRs",
                                       "Incomplete_IDs",
                                       "Incomplete_dates",
                                       "Incomplete_names",
                                       "Incomplete_material"),
                              Description = c("Total number of samples registered in Avenio_runs.xlsx",
                                              "Number of samples with complete info and able to be included in AVENIO_results_patients.rds",
                                              "Samples containing all the required information",
                                              "Number of samples with each designated type of material",
                                              "Number of samples with each designated sample type",
                                              "Runs entered in Avenio_runs.xlsx but is not present in AVENIO_results_patients.rds",
                                              "CPRs entered in Avenio_runs.xlsx but is not present in AVENIO_results_patients.rds",
                                              "Incomplete Run IDs",
                                              "Samples where dates are missing or wrongly formatted",
                                              "Samples where the sample name, project name or name in project is missing",
                                              "Samples where the material information is missing"))
    df <- df_current %>% 
        dplyr::left_join(explanation,by = "Name")
    return(df)
}