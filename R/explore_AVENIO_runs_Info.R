#' Highlights the parameters to be investigated from explore_AVENIO_runs()
#' 
#' This functions creates an overview of the different contents in the list
#' exported by explore_AVENIO_runs()
#' @importFrom dplyr `%>%` left_join
#' @return A `data.frame` with two variables containing the name and the
#' description of the entries in the explore_AVENIO_runs() output
#'  
#' @examples
#' explore_AVENIO_runs_Info()
#' @export
explore_AVENIO_runs_Info <- function(){
    current_format <- explore_AVENIO_runs(silent = TRUE)
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