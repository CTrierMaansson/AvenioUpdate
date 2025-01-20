#' Get a simple overview for a patient
#' 
#' This function takes the main list of NGS data for each patient
#' located at "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds"
#' and a CPR number (or other unique patient number used for the patient) and 
#' creates a simple overview of the mutations identified.
#' @importFrom dplyr left_join select mutate tibble `%>%` filter
#' @importFrom readxl read_xlsx
#' @importFrom BiocBaseUtils isScalarCharacter isScalarLogical
#' @param df_list The `list` of `data.frames` which contains all information on
#'  on the patients. This is read using 
#'  readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
#' @param  CPR `Character` string of the unique CPR number or other unique 
#'  number for the patient of interest. Available CPR numbers can be viewed
#'  using names(`df_list`)
#' @param synonymous `Boolean` determining whether synonymous mutations are kept
#'  in the output. If `TRUE` (default) synonymous mutations are kept in the 
#'  output. If `FALSE` the synonymous mutations are excluded. 
#' @return A `list` of length = 2. The first entry ("Patient CPR:") shows the 
#'  patient of interest CPR/ID number. The second entry ("Results") is a `tibble`
#'  with 14 variables. Each row corresponds to a mutation identified in any 
#'  sample. Important variables include sample_index which shows the project the
#'  patient is connected to, the name of the patient in that project and the 
#'  date of of the sample collection. Sample.ID is the name of the sample in the
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
#'  names(master)
#'  create_simple_output(df_list = master,
#'                       CPR = "1234567890",
#'                       synonymous = FALSE)
#' @export
create_simple_output <- function(df_list,CPR,synonymous = TRUE){
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
    if (!isScalarCharacter(CPR)) {
        stop("CPR has to be a character")
    }
    df <- df_list[[CPR]]
    if(is.null(df)){
        stop("CPR number is not in the dataset")
    }
    print("Reading AVENIO_runs.xlsx file")
    AVENIO_runs <- readxl::read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",
                             col_types = c(rep("guess",4),"date",rep("guess",6)))
    print("Selecting relevant sample")
    AVENIO_sele <- AVENIO_runs[AVENIO_runs$CPR == CPR,]
    AVENIO_sele <- AVENIO_sele[!is.na(AVENIO_sele$CPR),] %>% 
        dplyr::select(Sample_name, Sample_note, Material, Notes)
    print("Marging run and variant information")
    df <- df %>% 
        dplyr::left_join(AVENIO_sele, by = c("Sample.ID" = "Sample_name")) %>% 
        dplyr::select(sample_index,
               Mutation.Class,Gene,Amino.Acid.Change,
               Variant.Description,Flags,Allele.Fraction,
               Variant.Depth,Unique.Depth,Analysis.Name,Sample.ID,
               Sample_note, Material, Notes) %>%
        dplyr::mutate(Variant.Description = gsub(" variant","",Variant.Description)) %>% 
        dplyr::tibble()
    colnames(df) <- c("sample_index","Class","Gene","AA","Description",
                      "Flags","MAF","Variant_depth","Unique_depth","Analysis",
                      "Sample.ID","Sample_note","Material","Notes")
    if(synonymous){
        df <- df
    }
    else{
        print("Removing Synonymous variants")
        df <- df %>% 
            dplyr::filter(!grepl("Synonymous",Description))
    }
    return(list("Patient CPR:" = CPR,"Results" = df))
}