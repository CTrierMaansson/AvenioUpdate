#' @noRd
file_exist_func <- function(save_as_path,obj){
    `%ni%` <- Negate(`%in%`)
    f_exist <- file.exists(save_as_path)
    if(f_exist){
        m <- NA
        while(m %ni% c("1","2")){
            m <- readline(
                prompt = message(
                    paste0("The filename and path defined as\n",save_as_path,
                           " does exist.",
                           "\nAre you sure you want to override the file? \n", 
                           "1: Yes\n",
                           "2: No")))
            if(m %ni% c("1","2")){
                message("You have to write the options 1 or 2:")
            }
        }
        if(m == "1"){
            message(paste0("Saving results as:\n",
                           save_as_path))
            saveRDS(object = obj,
                    file = save_as_path)
        }
        if(m == "2")
            message("Select another filename using " ,
                    "save_as = <'file_name'>\n",
                    "Or do not save the file")
    }
    if(!f_exist){
        message(paste0("The filename ",save_as_path,"\n",
                       "Does not exist\n",
                       "\n",
                       "Using that path to save the file"))
        saveRDS(object = obj,
                file = save_as_path)
    }
}
