#' Highlights the parameters to be investigated from result_stats()
#' 
#' This functions creates an overview of the different contents in the list
#' exported by result_stats()
#' @importFrom dplyr `%>%` left_join
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/"
#' @return A `data.frame` with two variables containing the name and the
#' description of the entries in the result_stats() output
#'  
#' @examples
#' result_stats_Info()
#' @export
result_stats_Info <- function(synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    if (!dir.exists(synology_path)) {
        stop("The path entered as synology_path does not exist")
    }
    nchar_path <- nchar(synology_path)
    if(substr(synology_path,nchar_path,nchar_path) != "/"){
        stop("The synology_path has to end with a '/'")
    }
    current_format <- result_stats(silent = TRUE,
                                   synology_path = synology_path)
    df_current <- data.frame(Name = names(current_format))
    explanation <- data.frame(Name = c("Basestats",
                                       "Projectstats",
                                       "Missing",
                                       "All_mutations",
                                       "Relevant_SNV",
                                       "Relevant_INDEL",
                                       "BC_in_plasma",
                                       "Lengths",
                                       "Depths",
                                       "Reads",
                                       "On_target"),
                              Description = c("Number of samples, runs, patients, different materials etc.",
                                              "Number of patients and samples in each project",
                                              "Samples present in AVENIO_runs.xlsx but not present in the results data set",
                                              "Number of times each gene is mutated across all samples",
                                              "Number of SNVs detected not classified as BC or synonymous mutations",
                                              "Number of INDELs detected not classified as BC or synonymous mutations",
                                              "Number of times each gene is mutated in plasma but classified as BC mutation",
                                              "Distribution of fragmentlengths across different materials",
                                              "Distribution of unique depths across different materials",
                                              "Distribution of mapped reads across different materials",
                                              "Distribution of on target percents across different materials"))
    df <- df_current %>% 
        dplyr::left_join(explanation,by = "Name")
    return(df)
}
