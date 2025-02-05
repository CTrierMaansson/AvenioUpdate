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
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AajqaHzI9_hCH4wYCbDAIYx4")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABE9wiNEoilBeYhlV1UO3VRC")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AMw9wjlnV4lHRanHC8-jUnoh")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ANei_Hsi__lGmaJ1k_Ji_d3O")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AIH-1ryzAnFEMJ9R9-KI5_kA")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXP63fmerZZIn75npxDYZtYT")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ADx0jDAXEJ9Ho7oP9cWlxCB-")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKr3CLJPqtxPgYm8EWbJ8z_4")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AE5hBlinjEZN6ZZ1L7898Rzc")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXU04_NAKf5P3I1SrV_5x7oM")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AAMd0bIaeDJOcqthTfVBekrg")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AFzPNBUnw3NMOKplSiKRCqMU")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABeIKygVYa9ApLIqFBhxueub")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ACi6KICyf0tNSK8b3m4op5O0")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AAMz3RV2Tv1MrI1oA5zynUyG")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AHdh6Q0hlCRO263hSjLFJNQ7")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AR_CnOe8BaRPM71AHphnFT5J")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-Ado7dswg-kdHM5Fwb6dQYWOu")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKneBpDc1Q1Ln7Ev-5Iz53rK")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AGwQDzxcqrdDlpzlcS10s1En")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AU-GJVH_R-RAXa0UUil1MCik")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AYO0H3fNbPBEOKwJA9dv9IH0")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AMzRkv_uPpNMO6_XFokWfBuT")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AVhAZ6wwwchE4IIQa6lOQhN4")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ADNC-IlUdsROnLSz2AnnTpGz")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AGU232uEYZZKe4oI95IrJu7o")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AU-x-KCqRchJPaDwb-utcvvO")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AdA1JN9HI9ZND4LXP8XCfvRZ")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXVif6nxsixGLZDnsTvDtQLo")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= 30)
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ALsopAc5KC1Bx7pJLiHZ-eJR")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ASZMkYMxgytHgov13wWty29S")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AVhoQhsYxHRPKbWD3ZTc9uCa")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AEepGWqCaABI04g5DQESU97x")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ABSYoqHZHR9Gvo_Z0scoIQ8w")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AaRZMuj0NJtNnpZW8SO2dTy7")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AJD93ZaVEQxGcbVGcJdZZW94")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AdYl2rnCR1ZO5JLxgqP9fouS")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKorzNtokqVHNY5ZPcDPHeY3")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AY-MEcYGsnlEp70J19fZmwOJ")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AbZ6pmy7TndM4JOHXTcTCjde")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AeUIQZkyte1MjIyV0h0MS0aI")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AIjmFys4qs9EcrAcPI4IJxkK")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AaJsBOVJsSFO2onZbHBDtn0C")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKY-0douAMNBU7L0NimzcmCS")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ACHSPJ6rFXpCfKO7zOQZOWbz")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AHKsNmy-AFFNupR0Uot9fuks")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AZjnWu5wiDlJapI4CxvEuFnk")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AB2Pm9gToCJPGpSHjZEl7O5I")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AB3ubPmTO1hKvpB3SpZ0JSCj")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AceclBC2q4VHcKhBOYoHrQo9")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AWdzkUCGn0JHYp3dMNlyPyI0")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ASpCbtIIc89Mt4--IHz8OnPB")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-APzK4O25DCVD5YVyf63K1oGf")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>%
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AW7bfwpDfPRK5ogAOE44LCq_")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AGUZtNtuSiZEW4WOk9lYZtyN")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

# master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
# master_list <- AvenioUpdate::add_run_to_list(
#     master_list = master_list,
#     Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AfGU5cCRLV5CjaCoDbBjdzGN")
# analyses <- AvenioUpdate::included_analyses(master_list)
# print(analyses[[1]], n= nrow(analyses[[1]]))
# AvenioUpdate::result_stats()$Basestats
# tmp <- AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
#     filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AT5xbQo2WQtI2oi_dtrom7dM")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AFWAGp8DS6lNHrlt-I_azKUu")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ANEzevwHG-lNYLnNIR6aTGu9")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AIKiD4XtSL1NtbWtL5wIfq4m")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AXWJJwFGdYJBjbqclsNsl3UA")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ASP86HuqkG5G5KrN1nKqOGSh")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AZL9sfd4ALRFeqYfyNABLq04")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AQ53QFFOszBJZqx056aHcp4g")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKwjaEnnF1dAKKFcMt9igXXX")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AEHqc78-RkhKwqARpiG9V9Qi")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-Aai5AQdklPpMlJfseoQaBysU")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AauPU7MsWE1Dx51FLO3dJg33")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ARBybyYWlGJHUbhbQPibrXGI")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AEjxrRrIVjJGzbUFckTwpTpp")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AGB9eRbQCw9Bj63KZctA4DQe")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ARPd_ATt_oFBNrmQ7Yk_8YCo")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AYPCtGCABftK1I8Ej2WnbTc3")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-ARYBp-g_dE9Gra-ukrIXW-zc")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKoN6VzAOI9HgKcaCj1z9S8v")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AAOnQVP9mUBJk7Ij60bsjc6_")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AEusTCGgvJ1NrKQysOxTwSBe")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")


master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AKUURPIxRKxM5LcReW5cvFO9")
analyses <- AvenioUpdate::included_analyses(master_list)
print(analyses[[1]], n= nrow(analyses[[1]]))
AvenioUpdate::result_stats()$Basestats
AvenioUpdate::explore_AVENIO_runs(Info = "Unincluded_analyses") %>% 
    filter(Project != "SOLUCOM")
