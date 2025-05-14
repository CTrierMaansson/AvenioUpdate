#' Generates an overview different stats in the dataset
#' 
#' This function takes the AVENIO_results_patients.rds file and generates
#' different explorative stats on the dataset
#' @importFrom dplyr left_join select mutate `%>%` filter arrange count
#' @importFrom readxl read_xlsx
#' @importFrom BiocBaseUtils isScalarCharacter isScalarLogical
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
#' @param  Info `Character` string with the name of the information that should.
#' be extracted. Possible options can be viewed with result_stats_Info().
#' If `NULL` (default) all information is returned as a `list`.
#' @param silent `Boolean` determining information messages should be displayed 
#' If `FALSE` (default) messages are displayed. If `TRUE` the messages are
#' not printed. 
#' @param synology_path `Character` string with the full path to the directory 
#'  containing AVENIO_runs.xlsx & AVENIO_keys.rds. Default is 
#'  "//Synology_m1/Synology_folder/AVENIO/" 
#' @return A `list` with different stats on AVENIO_results_patients.rds
#' (`Info` = `NULL`) or the specific object as determined by `Info`.
#' 
#' @examples
#' results <- result_stats()
#' results$All_mutations
#' result_stats(Info = "Basestats")
#' @export
result_stats <- function(Info = NULL, silent = FALSE,
                         synology_path = "//Synology_m1/Synology_folder/AVENIO/"){
    `%ni%` <- Negate(`%in%`)
    if(!is.null(Info)){
        if (!isScalarCharacter(Info)) {
            stop("Info has to be a character or NULL")
        }
    }
    if (!isScalarLogical(silent)) {
        stop("silent has to be a TRUE or FALSE")
    }
    if (!dir.exists(synology_path)) {
        stop("The path entered as synology_path does not exist")
    }
    nchar_path <- nchar(synology_path)
    if(substr(synology_path,nchar_path,nchar_path) != "/"){
        stop("The synology_path has to end with a '/'")
    }
    if(!silent){
        message("Reading AVENIO_results_patients.rds")
    }
    master_list <- readRDS(paste0(synology_path,"AVENIO_results_patients.rds"))
    if(!silent){
        message("Creating combined data.frame")
    }
    combined_df <- do.call(rbind,master_list)
    if(!silent){
        message("Reading AVENIO_runs.xlsx")
    }
    rownames(combined_df) <- NULL
    Avenio_runs <- readxl::read_xlsx(
        paste0(synology_path,"AVENIO_runs.xlsx"),
        col_types = c(rep("guess",4),"date",rep("guess",6))) %>% 
        dplyr::mutate(date_check = lubridate::ymd(Sample_date)) %>% 
        dplyr::filter(CPR %in% names(master_list)) %>% 
        select(CPR,Name_in_project, Project, Sample_date,
               Run_name, Run_ID, Sample_name, 
               Sample_note, Material,date_check) %>% 
        na.omit()
    if(!silent){
        message("Creating sample_index variable")
    }
    Avenio_runs <- Avenio_runs %>% 
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
    if(!silent){
        message("Combining run and patient information")
    }
    combined_df <-  combined_df %>% 
        dplyr::left_join(Avenio_runs, by = "sample_index")
    gene_df <- combined_df %>% 
        dplyr::select(sample_index,Gene,Mutation.Class,Variant.Description,Flags,Material) %>% 
        unique()
    gene_df_complete <- gene_df %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    if(!silent){
        message("Extracting basic statistics")
    }
    basestats <- data.frame(stat = c("Patients",
                                     "Samples",
                                     "Runs",
                                     "Mutations"),
                            n = c(length(names(master_list)),
                                  length(unique(combined_df$sample_index)),
                                  length(unique(combined_df$Run_ID)),
                                  sum(gene_df_complete$n)))
    mat_df <- combined_df %>% 
        dplyr::select(sample_index,Material) %>% 
        dplyr::group_by(Material) %>% 
        unique() %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    note_df <- combined_df %>% 
        dplyr::select(sample_index,Sample_note) %>% 
        dplyr::group_by(Sample_note) %>% 
        unique() %>% 
        dplyr::count() %>% 
        dplyr::filter(Sample_note != "BC") %>% 
        dplyr::arrange(desc(n))
    colnames(mat_df)[1] <- "stat"
    colnames(note_df)[1] <- "stat"
    basestats <- rbind(basestats, mat_df, note_df)
    
    if(!silent){
        message("Finding unincluded samples")
    }
    total_df <- do.call(rbind,master_list)
    rownames(total_df) <- NULL
    added_samples <- unique(total_df$sample_index)
    xlsx_samples <- unique(Avenio_runs$sample_index)
    missing_samples <- xlsx_samples[xlsx_samples %ni% added_samples]
    problematic_df <- Avenio_runs %>% 
        dplyr::filter(sample_index %in% missing_samples) %>% 
        dplyr::select(sample_index, Run_ID,Sample_name) %>% 
        as.data.frame()
    if(!silent){
        message("Extracting project statistics")
    }
    project_df_patients <- combined_df %>% 
        dplyr::select(CPR,Project) %>% 
        unique() %>% 
        dplyr::group_by(Project) %>% 
        dplyr::count()
    colnames(project_df_patients)[2] <- "n_patients"
    project_df_samples <- combined_df %>% 
        dplyr::select(sample_index,Project) %>% 
        unique() %>% 
        dplyr::group_by(Project) %>% 
        dplyr::count()
    colnames(project_df_samples)[2] <- "n_samples"
    project_df <- project_df_patients %>% 
        dplyr::left_join(project_df_samples, by = "Project") %>% 
        dplyr::ungroup()
    if(!silent){
        message("Extracting gene mutation statistics")
    }
    gene_df <- combined_df %>% 
        dplyr::select(sample_index,Gene,Mutation.Class,Variant.Description,Flags,Material) %>% 
        unique()
    gene_df_complete <- gene_df %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    gene_df_SNV_no_BC <- gene_df %>% 
        dplyr::filter(Material != "BC") %>% 
        dplyr::filter(!grepl("Synonymous",Variant.Description)) %>% 
        dplyr::filter(!grepl("BC",Flags)) %>%
        dplyr::filter(Mutation.Class == "SNV") %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    gene_df_INDEL_no_BC <- gene_df %>% 
        dplyr::filter(Material != "BC") %>% 
        dplyr::filter(!grepl("Synonymous",Variant.Description)) %>% 
        dplyr::filter(!grepl("BC",Flags)) %>%
        dplyr::filter(Mutation.Class == "INDEL") %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    BC_in_cfDNA <- gene_df %>% 
        dplyr::filter(Material != "BC") %>% 
        dplyr::filter(grepl("BC",Flags)) %>%
        dplyr::filter(Mutation.Class %in% c("INDEL","SNV")) %>% 
        dplyr::group_by(Gene) %>% 
        dplyr::count() %>% 
        dplyr::arrange(desc(n))
    if(!silent){
        message("Extracting NGS base statistics")
    }
    size_depth_df <- combined_df %>% 
        dplyr::select(Unique.Depth.Median,Mapped.Reads,
                      On.Target.Rate,Fragment.Length.Median,sample_index,Material) %>% 
        unique() %>% 
        dplyr::mutate(Unique.Depth.Median = as.numeric(Unique.Depth.Median)) %>% 
        dplyr::mutate(Mapped.Reads = as.numeric(Mapped.Reads)) %>% 
        dplyr::mutate(On.Target.Rate = gsub("%","",On.Target.Rate)) %>% 
        dplyr::mutate(On.Target.Rate = as.numeric(On.Target.Rate)) %>% 
        dplyr::mutate(Fragment.Length.Median = as.numeric(Fragment.Length.Median)) %>% 
        dplyr::select(-sample_index)
    for(i in 1:length(unique(size_depth_df$Material))){
        mat <- unique(size_depth_df$Material)[i]
        size_depth_df_mat <- size_depth_df %>% 
            dplyr::filter(Material == mat)
        hist_reads <- hist(size_depth_df_mat$Mapped.Reads,breaks = 20,plot = F)
        hist_target <- hist(size_depth_df_mat$On.Target.Rate,breaks = 20,plot = F)
        hist_length <- hist(size_depth_df_mat$Fragment.Length.Median,breaks = 20,plot = F)
        hist_depth <- hist(size_depth_df_mat$Unique.Depth.Median,breaks = 20,plot = F)
        df_target_mat <- data.frame(Target = hist_target$mids,
                                    n = hist_target$counts,
                                    Material = mat)
        df_reads_mat <- data.frame(Reads = hist_reads$mids,
                                    n = hist_reads$counts,
                                    Material = mat)
        df_length_mat <- data.frame(Length = hist_length$mids,
                                    n = hist_length$counts,
                                    Material = mat)
        df_depth_mat <- data.frame(Depth = hist_depth$mids,
                                    n = hist_depth$counts,
                                    Material = mat)
        if(i == 1){
            df_target <- df_target_mat
            df_reads <- df_reads_mat
            df_length <- df_length_mat
            df_depth <- df_depth_mat
        }else{
            df_target <- rbind(df_target,df_target_mat)
            df_length <- rbind(df_length,df_length_mat)
            df_reads <- rbind(df_reads,df_reads_mat)
            df_depth <- rbind(df_depth,df_depth_mat)
        }
    }
    if(!silent){
        message("Formatting output")
    }
    res <- list(Basestats = basestats,
                Projectstats = project_df,
                Missing = problematic_df,
                All_mutations = gene_df_complete,
                Relevant_SNV = gene_df_SNV_no_BC,
                Relevant_INDEL = gene_df_INDEL_no_BC,
                BC_in_plasma = BC_in_cfDNA,
                Lengths = df_length,
                Depths = df_depth,
                Reads = df_reads,
                On_target = df_target)
    if(!is.null(Info)){
        if(Info %ni% names(res)){
            stop("Info has to a a viable name of the output list. Relevant entries can be viewed with result_stats_Info()")
        }
        if(!silent){
            message(paste0("Outputting: ",Info))
        }
        output <- res[[Info]]
    }
    else{
        output <- res
    }
    return(output)
}