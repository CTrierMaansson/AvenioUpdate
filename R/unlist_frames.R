#' Create one big data.frame from list of patient data.frames
#' 
#' This function takes the main list of data.frames containing the mutation 
#' information and creates a single data.frame with the same information, 
#' where each individual is concatenated
#' @importFrom dplyr `%>%` left_join bind_rows
#' @param master_list The `list` of `data.frames` which contains all information on
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @return A `data.frame` with the all samples
#' @examples
#' master <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' master_df <- unlist_frames(master_list = master)
#' print(master_df)
#' @export
unlist_frames <- function(master_list){
    if (!is.list(master_list)) {
        stop("df_list has to be a list")
    }
    for(i in 1:length(master_list)){
        if(!is.data.frame(master_list[[i]])){
            stop(paste0("Entry ",i," in the master_list is not a data.frame"))
        }
    }
    res <- do.call(dplyr::bind_rows,master_list) %>% 
        unique()
    sample_nam <- do.call(dplyr::bind_rows, 
                      lapply(names(master_list), 
                             function(nm) {
        data.frame(CPR = nm, 
                   sample_index = master_list[[nm]][[1]])
                                 })) %>% 
        unique()
    res <- res %>% 
        left_join(sample_nam, by = "sample_index")
    return(res)
}


