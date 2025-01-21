#If the "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds"
#file should get corrupted or otherwise get lost/missing this script
#can recreate the file which has been updated as of 21-01-2025

#IMPORTANT! DO NOT RUN THIS UNLESS ABSOLUTELY NECESSARY!!!

zip_path1 <- "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AbpTmAnr_h9DPL1skMvyje5U"
library(readxl)
library(dplyr)
library(stringr)
`%ni%` <- Negate(`%in%`)

add_samples_mod <- function(Directory){
    Run_ID <- gsub("//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-","",Directory)
    zip_file <- paste0(Directory,"/","all_CSVs-",Run_ID,".zip")
    sample_file <- paste0(Directory,"/SampleMetrics.csv")
    run_ID_short <- substr(Run_ID,1,8)
    unzip(zipfile = zip_file,
          exdir = Directory)
    Filter_variants <- paste0(Directory,"/","VariantMetrics-UserFilterSet.csv")
    variants <- read.csv(Filter_variants)
    AVENIO_runs <- read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",col_types = c(rep("guess",4),"date",rep("guess",6)))
    AVENIO_runs_select <- AVENIO_runs %>% 
        filter(grepl(run_ID_short,Run_ID))
    variants_select <- variants %>% 
        filter(Sample.ID %in% AVENIO_runs_select$Sample_name)
    samples_with_variants <- unique(variants_select$Sample.ID)
    samples_of_interest <- unique(AVENIO_runs_select$Sample_name)
    if(length(samples_with_variants) != length(samples_of_interest)){
        sample_metrics <- read.csv(sample_file) %>% 
            filter(Sample.ID %in% samples_of_interest) %>% 
            filter(Sample.ID %ni% samples_with_variants)
        sample_metrics <- sample_metrics %>% 
            select(colnames(sample_metrics)[colnames(sample_metrics) %in% colnames(variants_select)])
        sele_cols <- colnames(variants_select)[colnames(variants_select) %in% colnames(sample_metrics)]
        for(i in 1:length(sele_cols)){
            if(class(sample_metrics[[sele_cols[i]]]) != class(variants_select[[sele_cols[i]]])){
                classes <- c(class(sample_metrics[[sele_cols[i]]]),class(variants_select[[sele_cols[i]]]))
                if("character" %in% classes){
                    sample_metrics[[sele_cols[i]]] <- as.character(sample_metrics[[sele_cols[i]]])
                    variants_select[[sele_cols[i]]] <- as.character(variants_select[[sele_cols[i]]])
                }
                if("numeric" %in% classes){
                    sample_metrics[[sele_cols[i]]] <- as.numeric(sample_metrics[[sele_cols[i]]])
                    variants_select[[sele_cols[i]]] <- as.numeric(variants_select[[sele_cols[i]]])
                }
            }
        }
        combined_samples <- full_join(variants_select,sample_metrics,copy = TRUE)
    }
    else{
        combined_samples <- variants_select
    }
    if(!"sample_index" %in% colnames(AVENIO_runs_select)){
        AVENIO_runs_select <- AVENIO_runs_select %>% 
            mutate(sample_index = paste0(Project,"_",
                                         Name_in_project,"_",
                                         substr(str_split_i(as.character(Sample_date),"-",1),3,4),
                                         str_split_i(as.character(Sample_date),"-",2),
                                         str_split_i(as.character(Sample_date),"-",3))) %>% 
            mutate(sample_index = ifelse(Material != "cfDNA",paste0(sample_index,"_",Material),sample_index))
    }
    AVENIO_runs_select_merge <- AVENIO_runs_select %>% 
        dplyr::select(Sample_name,sample_index)
    combined_samples <- combined_samples %>% 
        left_join(AVENIO_runs_select_merge,by = c("Sample.ID"="Sample_name")) %>% 
        relocate(sample_index,.before = Sample.ID) %>% 
        arrange(sample_index,Gene)
    return(combined_samples)
}
combined_samples1 <- add_samples_mod(zip_path1)

create_df_list_mod <- function(sample_df){
    samples <- unique(sample_df$Sample.ID)
    AVENIO_runs <- read_xlsx("//Synology_m1/Synology_folder/AVENIO/AVENIO_runs.xlsx",col_types = c(rep("guess",4),"date",rep("guess",6)))
    AVENIO_runs_select <- AVENIO_runs %>% 
        filter(Sample_name %in% samples)
    unique_patients <- unique(AVENIO_runs_select$CPR)
    df_list <- list()
    for (i in 1:length(unique_patients)){
        CPR_patient <- unique_patients[i]
        patient_df <- AVENIO_runs_select %>% 
            filter(CPR == CPR_patient)
        patient_variants <- sample_df %>% 
            filter(Sample.ID %in% patient_df$Sample_name)
        df_list[[i]] <- patient_variants 
        names(df_list)[i] <- CPR_patient
    }
    return(df_list)
}
master_list <- create_df_list_mod(combined_samples1)

master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABE9wiNEoilBeYhlV1UO3VRC")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AMw9wjlnV4lHRanHC8-jUnoh")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ANei_Hsi__lGmaJ1k_Ji_d3O")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AIH-1ryzAnFEMJ9R9-KI5_kA")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXP63fmerZZIn75npxDYZtYT")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ADx0jDAXEJ9Ho7oP9cWlxCB-")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKr3CLJPqtxPgYm8EWbJ8z_4")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AE5hBlinjEZN6ZZ1L7898Rzc")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXU04_NAKf5P3I1SrV_5x7oM")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AAMd0bIaeDJOcqthTfVBekrg")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AFzPNBUnw3NMOKplSiKRCqMU")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABeIKygVYa9ApLIqFBhxueub")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ACi6KICyf0tNSK8b3m4op5O0")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AAMz3RV2Tv1MrI1oA5zynUyG")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AHdh6Q0hlCRO263hSjLFJNQ7")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AR_CnOe8BaRPM71AHphnFT5J")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-Ado7dswg-kdHM5Fwb6dQYWOu")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
names(master_list)

