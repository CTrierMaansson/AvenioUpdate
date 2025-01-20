#' Get an overview of the analyses included in the dataset
#' 
#' This function takes the main dataset and creates an overview of how many 
#' samples from each NGS run is included in the dataset
#' @importFrom dplyr select `%>%` filter group_by left_join count
#' @param master_list The `list` of `data.frames` which contains all information
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @return A `list` of length = 2. The first entry ("Overview:") 
#' is a `tibble` with three variables giving an overview of the number of 
#' patients (n) from each NGS run included in the dataset. The second entry is a
#' `list` of `data.frames` where each `data.frame` represents the samples from 
#' the specific NGS run. Each `data.frame` is named according to the Analysis.ID
#'  
#' @examples
#' master <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' analyses <- included_analyses(master)
#' analyses[[1]]
#' names(analyses[[2]])
#' analyses[[2]][["AbpTmAnr_h9DPL1skMvyje5U"]]
#' @export
included_analyses <- function(master_list){
    if (!is.list(master_list)) {
        stop("df_list has to be a list")
    }
    for(i in 1:length(master_list)){
        if(!is.data.frame(master_list[[i]])){
            stop(paste0("Entry ",i," in the master_list is not a data.frame"))
        }
    }
    combined_df <- do.call(rbind,master_list) %>% 
        dplyr::select(sample_index,Sample.ID,Analysis.Name,Analysis.ID) %>% 
        unique()
    rownames(combined_df) <- NULL
    unique_IDs <- unique(combined_df$Analysis.ID)
    res_list <- list()
    for(i in 1:length(unique_IDs)){
        ID <- unique_IDs[i]
        df_ID <- combined_df %>% 
            dplyr::filter(Analysis.ID == ID)
        res_list[[i]] <- df_ID
        names(res_list)[i] <- ID
    }
    combined_sele <- combined_df %>% 
        dplyr::select(Analysis.ID, Analysis.Name) %>% 
        unique()
    combined_sum <- combined_df %>% 
        dplyr::group_by(Analysis.ID) %>% 
        dplyr::count() %>% 
        dplyr::left_join(combined_sele, by = "Analysis.ID") %>% 
        ungroup()
    return(list("Overview:" = combined_sum,
                "Details:" = res_list))
}
included_analyses(master_list)
