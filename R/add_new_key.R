#' Add a new key to Project, Sample_note, or Material
#' 
#' This function can add new possible keys for Project, Sample_note, or Material
#' in the AVENIO_runs.xlsx file. 
#' @importFrom BiocBaseUtils isScalarCharacter
#' @importFrom readxl read_xlsx
#' @param key `Character` of the key name you want to add as possible entry
#' @param  variable `Character` of the variable name in the AVENIO_runs.xlsx
#' file you want to add the key to.
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/"  
#' @return Nothing is returned. The key is just added to the list of possible
#' entries
#'  
#' @examples
#' add_new_key(key = "test",
#'             variable = "Project")
#' @export
add_new_key <- function(key,variable,
                        synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    if (!isScalarCharacter(key)) {
        stop("key has to be a character")
    }
    if (!isScalarCharacter(variable)) {
        stop("variable has to be a character")
    }
    `%ni%` <- Negate(`%in%`)
    Avenio_runs <- readxl::read_xlsx(
        paste0(synology_path,"AVENIO_runs.xlsx"),
        col_types = c(rep("guess",4),"date",rep("guess",6)))
    keys <- readRDS(paste0(synology_path,"AVENIO_keys.rds"))
    if(variable %ni% names(keys)){
        key_names <- paste(names(keys),collapse = ", ")
        stop(paste0("variable name: '", variable,
                    "' is invalid, possible entries include:\n",
                    key_names))
    }
    key_entries <- keys[[variable]]
    if(key %in% key_entries){
        stop(paste0("key name: '", key,
                    "' already exists for the variable: '",
                    variable,"' and will not be added again.\nTerminating"))
    }
    if(variable == "Project"){
        if(grepl("_",key)){
            stop("The Project variable cannot contain entries with '_'.\n",
                 "Please change the name of the key: '",key,"'.")
        }
    }
    key_entries <- c(key_entries,key)
    keys_new <- keys
    keys_new[[variable]] <- key_entries
    tests <- c()
    for(i in 1:length(names(keys))){
        var <- names(keys)[i]
        ent_old <- keys[[var]]
        ent_new <- keys_new[[var]]
        tests[i] <- length(ent_new) >= length(ent_old)
    }
    if(all(tests)){
        message(paste0("Adding the new key: '",key,"' for the variable: '",
                       variable,"'."))
        saveRDS(keys_new,paste0(synology_path,"AVENIO_keys.rds"))
        message("DONE!")
    }
}