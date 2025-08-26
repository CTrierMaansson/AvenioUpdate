#' @noRd
#' @importFrom dplyr mutate `%>%` filter
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
create_sample_index <- function(df){
    df <- df %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date)) %>% 
        dplyr::filter(!is.na(date_check)) %>%
        dplyr::mutate(sample_index = paste0(Project,"_",
                                            Name_in_project,"_",
                                            substr(stringr::str_split_i(
                                                as.character(Sample_date),"-",1),
                                                3,4),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",2),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",3))) %>% 
        dplyr::mutate(
            sample_index = ifelse(Material == "cfDNA",
                                  sample_index,
                                  ifelse(Material == "tissue",
                                         paste0(sample_index,"_",
                                                Material,"_",
                                                Sample_note),
                                         ifelse(Material == "reanalyze",
                                                ifelse(Sample_note == "BC",
                                                       paste0(sample_index,"_",
                                                              Sample_note,"_",
                                                              Material),
                                                       paste0(sample_index,"_",
                                                              Material)),
                                                paste0(sample_index,"_",
                                                       Material)))))
    return(df)
}