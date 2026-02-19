#' @noRd
save_as_func <- function(save_as_path){
    `%ni%` <- Negate(`%in%`)
    n <- NA
    while(n %ni% c("1","2")){
        n <- readline(
            prompt = message(
                paste0("Select if output should be saved:\n",
                       "1: Save object\n",
                       "2: Do not save object")))
        if(n %ni% c("1","2")){
            message("You have to write the options 1 or 2:")
        }
    }
    return(n)
}
