#If the "//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds"
#file should get corrupted or otherwise get lost/missing this script
#can recreate the file which has been updated as of 21-01-2025

#IMPORTANT! DO NOT RUN THIS UNLESS ABSOLUTELY NECESSARY!!!

master_list <- list()

master_list <- AvenioUpdate::add_run_to_list(
    master_list = master_list,
    Directory = "//Synology_m1/Synology_folder/AVENIO/AVENIO_results/Plasma-AbpTmAnr_h9DPL1skMvyje5U")
analyses <- AvenioUpdate::included_analyses(master_list)
analyses[[1]]

master_list <- readRDS("//Synology_m1/Synology_folder/AVENIO/AVENIO_results_patients.rds")
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

