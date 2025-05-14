#' @noRd
#' @importFrom dplyr filter `%>%` mutate select bind_rows left_join arrange
#' @importFrom stringr str_split_i
#' @importFrom lubridate ymd
#' @importFrom DNAfusion EML4_ALK_analysis find_variants
add_DNAfusion_res <- function(df,
                              repo_path,
                              sample_info){
    unique_samples <- unique(df$sample_index)
    sample_info <- sample_info %>% 
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
            sample_index = ifelse(Material != "cfDNA",
                                  paste0(sample_index,"_",Material),
                                  sample_index))
    for(i in 1:length(unique_samples)){
        sample_ind <- unique_samples[i]
        sample_info <- df %>% 
            dplyr::filter(sample_index == sample_ind)
        sampleID <- unique(sample_info$Sample.ID)
        panelID <- unique(sample_info$Panel)
        panelID <- ifelse(panelID == "Surveillance v2","Surveillance-v2",panelID)
        bam_file <- paste0(repo_path,"/",sampleID,"-",panelID,"/PosDeduped-",sampleID,".bam")
        DNAfusion_res <- DNAfusion::EML4_ALK_analysis(bam_file)
        if(length(DNAfusion_res) != 0){
            clipped_reads <- max(c(as.numeric(max(DNAfusion_res$breakpoint_ALK)),
                                   as.numeric(max(DNAfusion_res$breakpoint_EML4))),
                                 na.rm = T)
            if(is.character(names(which.max(DNAfusion_res$breakpoint_ALK)))){
                logic_ALK_breakpoint = T
                fusion_position_1 <- paste0("chr2:",
                                            names(which.max(DNAfusion_res$breakpoint_ALK)))
            }
            else{
                fusion_position_1 <- paste0("chr2:","NA")
                logic_ALK_breakpoint = F
            }
            if(is.character(names(which.max(DNAfusion_res$breakpoint_EML4)))){
                logic_EML4_breakpoint = T
                fusion_position_2 <- paste0("chr2:",
                                            names(which.max(DNAfusion_res$breakpoint_EML4)))
            }
            else{
                fusion_position_2 <- paste0("chr2:","NA")
                logic_EML4_breakpoint = F
            }
            fusion_position <- paste0(fusion_position_1,";",fusion_position_2)
            if(logic_ALK_breakpoint){
                MAF_num <- round((clipped_reads/DNAfusion_res$read_depth_ALK)*100,digits = 2)
                MAF <- paste0(MAF_num,"%")
                ALK_depth <- as.character(DNAfusion_res$read_depth_ALK)
            }
            else{
                MAF <- NA
                ALK_depth <- NA
            }
            if(logic_EML4_breakpoint){
                message(paste0("Establishing the EML4-ALK variant"))
                fusion_variant_res <- DNAfusion::find_variants(bam_file)
                if(is.list(fusion_variant_res)){
                    if(is.data.frame(fusion_variant_res)){
                        fusion_variant <- paste0(fusion_variant_res$Variant)
                        exon <- paste0("ALK_i",fusion_variant_res$Intron_ALK,
                                       ";EML4_i",fusion_variant_res$Intron_EML4)
                    }else{
                        fusion_variant <- "Not classified"
                        exon <- paste0("ALK_i",fusion_variant_res$ALK_intron,
                                       ";EML4_i",fusion_variant_res$EML4_intron)
                    }
                } 
                else{
                    fusion_variant <- fusion_variant_res
                    exon <- paste0("ALK_i",NA,
                                   ";EML4_i",NA)
                }
            }else{
                message(paste0("Could not establish EML4 breakpoint and \n",
                               "thefore cannot establish EML4-ALK variant"))
                exon <- paste0("ALK_i",NA,
                               ";EML4_i",NA)
                fusion_variant <- "Not classified"
            }
            res_df <- data.frame(sample_index = sample_ind,
                                 Sample.ID = sampleID,
                                 Flags = "DNAfusion",
                                 Mutation.Class = "FUSION",
                                 Gene = "ALK;EML4",
                                 Variant.Description = fusion_variant,
                                 Genomic.Position = fusion_position,
                                 Fusion.Partner.Gene.1 = "ALK",
                                 Fusion.Partner.Gene.2 = "EML4",
                                 Fusion.Breakpoint.1 = fusion_position_1,
                                 Fusion.Breakpoint.2 = fusion_position_2,
                                 Paired.Reads.Spanning.Fusion.Breakpoint = as.character(clipped_reads),
                                 Allele.Fraction = MAF,
                                 Variant.Depth = as.character(clipped_reads),
                                 Unique.Depth = ALK_depth,
                                 Exon.Number = exon)
            sample_info_unique <- sample_info %>% 
                dplyr::select(colnames(df)[c(1,4:15,47:84)]) %>% 
                unique() %>% 
                dplyr::left_join(res_df, by = "sample_index")
            sample_res <- sample_info %>% 
                dplyr::mutate(Variant.Depth = as.character(Variant.Depth)) %>% 
                dplyr::mutate(
                    Paired.Reads.Spanning.Fusion.Breakpoint = 
                        as.character(Paired.Reads.Spanning.Fusion.Breakpoint)) %>% 
                dplyr::mutate(Unique.Depth = as.character(Unique.Depth)) %>% 
                dplyr::bind_rows(sample_info_unique) %>% 
                dplyr::arrange(Gene) %>% 
                dplyr::filter(!is.na(Gene))
            df <- df %>% 
                dplyr::filter(sample_index != sample_ind) %>% 
                rbind(sample_res)
        }
    }
    df <- df %>% 
        dplyr::arrange(sample_index,Gene)
    return(df)
}
add_DNAfusion_res(df = samples,
                  repo_path = Directory,
                  sample_info = AVENIO_runs_select)
clipped_reads <- max(c(as.numeric(max(res$breakpoint_ALK)),
                       as.numeric(max(res$breakpoint_EML4))),
                     na.rm = T)
if(is.character(names(which.max(res$breakpoint_ALK)))){
    logic_ALK_breakpoint = T
    fusion_position_1 <- paste0("chr2:",
                                names(which.max(res$breakpoint_ALK)))
}else{
    fusion_position_1 <- paste0("chr2:","NA")
    logic_ALK_breakpoint = F
}
if(is.character(names(which.max(res$breakpoint_EML4)))){
    logic_EML4_breakpoint = T
    fusion_position_2 <- paste0("chr2:",
                                names(which.max(res$breakpoint_EML4)))
}else{
    fusion_position_2 <- paste0("chr2:","NA")
    logic_EML4_breakpoint = F
}
fusion_position <- paste0(fusion_position_1,";",fusion_position_2)
if(logic_ALK_breakpoint){
    MAF_num <- round((clipped_reads/res$read_depth_ALK)*100,digits = 2)
    MAF <- paste0(MAF_num,"%")
    ALK_depth <- as.character(res$read_depth_ALK)
}else{
    MAF <- NA
    ALK_depth <- NA
}
