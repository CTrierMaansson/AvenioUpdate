#' Highlights the parameters to be investigated from result_stats()
#' 
#' This functions creates an overview of the different contents in the list
#' exported by result_stats()
#' @importFrom dplyr `%>%` left_join
#' @return A `data.frame` with two variables containing the name and the
#' description of the entries in the result_stats() output
#'  
#' @examples
#' result_stats_Info()
#' @export
result_stats_Info <- function(){
    current_format <- result_stats(silent = TRUE)
    df_current <- data.frame(Name = names(current_format))
    explanation <- data.frame(Name = c("Basestats",
                                       "Projectstats"
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