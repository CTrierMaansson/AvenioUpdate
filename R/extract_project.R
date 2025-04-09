#' Get a simple overview for an entire project
#' 
#' This function takes the main list of NGS data for each patient
#' located at "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds"
#' and a project ID and creates a simple overview of the mutations identified
#' for the patients in that project
#' @importFrom dplyr left_join select mutate tibble `%>%` filter arrange
#' @importFrom readxl read_xlsx
#' @importFrom BiocBaseUtils isScalarCharacter isScalarLogical
#' @importFrom stringr str_split_i
#' @param df_list The `list` of `data.frames` which contains all information on
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @param  project `Character` string of the project name as it appears in
#' AVENIO_runs.xlsx. Accepted project names are available with 
#' `result_stats(Info = "Projectstats")` 
#' @param synonymous `Boolean` determining whether synonymous mutations are kept
#'  in the output. If `TRUE` (default) synonymous mutations are kept in the 
#'  output. If `FALSE` the synonymous mutations are excluded. 
#' @param simple `Boolean`. If `TRUE` the output is in a simple format,
#' similar to the format created by `create_simple_output()`. If `FALSE` 
#' (default) all information from the AVENIO .csv files are exported. 
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/" 
#' @return A `tibble`
#'  with 15 variables (`simple = TRUE`) or many variables (`simple = FALSE`).
#'  Each row corresponds to a mutation identified in any 
#'  sample. Important variables include sample_index which shows the project the
#'  patient is connected to, the name of the patient in that project and the 
#'  date of of the sample collection. CPR which is the CPR number of the patient,
#'  Sample.ID is the name of the sample in the
#'  AVENIO system. Sample_note is the note specifying what type of blood 
#'  sample the sample is in relation to treatment (BL, Tx, BC, or Unknown). 
#'  Material specifies what type of material has been analyzed (cfDNA, BC, 
#'  cfChIP, size-selection, etc.) The Flags variable shows if the mutation only
#'  has been found in the BAM file and not in the AVENIO output, and/or if the
#'  mutation also has been found in a BC sample if such sample has been run for
#'  the patient. If a sample has been analyzed but no mutation was identified 
#'  `NA` is in all relevant variables.
#'  
#' @examples
#'  master <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#'  result_stats(Info = "Projectstats")
#'  extract_project(df_list = master,
#'                  project = "Pembrolizumab",
#'                  synonymous = FALSE,
#'                  simple = TRUE)
#' @export
extract_project <- function(df_list,
                            project,
                            synonymous = TRUE,
                            simple = FALSE,
                            synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    `%ni%` <- Negate(`%in%`)
    if (!is.list(df_list)) {
        stop("df_list has to be a list")
    }
    for(i in 1:length(df_list)){
        if(!is.data.frame(df_list[[i]])){
            stop(paste0("Entry ",i," in the df_list is not a data.frame"))
        }
    }
    if (!isScalarLogical(synonymous)) {
        stop("synonymous has to be a TRUE or FALSE")
    }
    if (!isScalarLogical(simple)) {
        stop("simple has to be a TRUE or FALSE")
    }
    if (!isScalarCharacter(project)) {
        stop("project has to be a character")
    }
    if (!dir.exists(synology_path)) {
        stop("The path entered as synology_path does not exist")
    }
    nchar_path <- nchar(synology_path)
    if(substr(synology_path,nchar_path,nchar_path) != "/"){
        stop("The synology_path has to end with a '/'")
    }
    message("Reading AVENIO_runs.xlsx file")
    AVENIO_runs <- readxl::read_xlsx(paste0(synology_path,"AVENIO_runs.xlsx"),
                                     col_types = c(rep("guess",4),"date",rep("guess",6))) 
    AVENIO_runs <- AVENIO_runs %>% 
        dplyr::mutate(sample_index = paste0(Project,"_",
                                            Name_in_project,"_",
                                            substr(stringr::str_split_i(
                                                as.character(Sample_date),"-",1),
                                                3,4),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",2),
                                            stringr::str_split_i(
                                                as.character(Sample_date),"-",3)
        )
        ) %>% 
        dplyr::mutate(
            sample_index = ifelse(Material != "cfDNA",
                                  paste0(sample_index,"_",Material),
                                  sample_index))
    if(project %ni% AVENIO_runs$Project){
        stop("project is invalid. Does not appear in AVENIO_runs.xlsx\nSee available projects with result_stats(Info = 'Projectstats')")
    }
    message("Selecting relevant sample")
    AVENIO_runs_sele <- AVENIO_runs %>% 
        dplyr::filter(Project == project) %>% 
        dplyr::select(sample_index,CPR,Sample_note, Material, Notes) %>% 
        dplyr::filter(!is.na(sample_index))
    if(is.null(AVENIO_runs_sele)){
        stop("No patients are available in the selected project")
    }
    combined_df <- do.call(rbind,df_list)
    rownames(combined_df) <- NULL
    df <- combined_df %>% 
        unique() %>% 
        dplyr::filter(sample_index %in% AVENIO_runs_sele$sample_index) %>% 
        dplyr::mutate(Variant.Description = gsub(" variant","",Variant.Description)) %>% 
        as.data.frame() %>% 
        dplyr::left_join(AVENIO_runs_sele,by = "sample_index") %>% 
        filter(grepl("pt01_",sample_index))
    if(simple){
        message("Creating output in simple format")
        df <- df %>% 
            dplyr::select(sample_index, CPR,
                          Mutation.Class,Gene,Amino.Acid.Change,
                          Variant.Description,Flags,Allele.Fraction,
                          Variant.Depth,Unique.Depth,Analysis.Name,Sample.ID,
                          Sample_note, Material, Notes)
        colnames(df) <- c("sample_index","CPR","Class","Gene","AA","Description",
                          "Flags","MAF","Variant_depth","Unique_depth","Analysis",
                          "Sample.ID","Sample_note","Material","Notes")
        df <- df %>% 
            dplyr::arrange(sample_index,Gene)
        if(synonymous){
            df <- df
        }
        else{
            message("Removing Synonymous variants")
            df <- df %>% 
                dplyr::filter(!grepl("Synonymous",Description))
        }
    }
    else{
        df <- df %>% 
            dplyr::arrange(sample_index,Gene)
        if(synonymous){
            df <- df
        }
        else{
            message("Removing Synonymous variants")
            df <- df %>% 
                dplyr::filter(!grepl("Synonymous",Variant.Description))
        }
    }
    
    return(df)
}