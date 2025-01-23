#' @noRd
#' @importFrom dplyr tibble `%>%` filter
#' @importFrom stringr str_glue
#' @importFrom vcfR read.vcfR
wrapper_tjek_bam_snv<-function(df){
    
    df_varianter_fundet_i_bam<-dplyr::tibble(
        sample_1=c(0),sample_2=c(0),
        genomic_position=c(0),
        variant_depth=c(0),
        unique_depth=c(0))
    df_nested<-lav_nested_df_SNV(df)
    path_AVENIO_results<-"//Synology_m1/Synology_folder/AVENIO/AVENIO_results"
    
    for(i in 1:nrow(df_nested)){
        #get sample1 name:
        sample1<-df_nested$Sample[i]
        print(paste0("Processing ", sample1))
        
        #path til Bam fil sample i:
        analyseID<-df_nested$Analysis.ID[i]
        sampleID<-df_nested$Sample.ID[i]
        panelID<-df_nested$Panel[i]
        panelID<-ifelse(panelID == "Surveillance v2","Surveillance-v2",panelID)
        path_sample_mappe<-paste0(path_AVENIO_results,
                                  "/Plasma-",
                                  analyseID,
                                  "/",
                                  sampleID,
                                  "-",
                                  panelID)
        #Andre prøver fra samme patient:
        df_patient_other_samples<-df_nested%>%
            dplyr::filter(Patient==df_nested$Patient[i]&Sample!=sample1)
        
        #hvis der er kørt flere prøver på patienten, så tjekkes disse en efter en:
        if (nrow(df_patient_other_samples) == 0) 
        {print(paste0("No other samples for the patient than ", sample1))}
        else{
            #loop over andre prøver:
            for (j in 1:nrow(df_patient_other_samples)){
                #get sample2 name:
                sample2<-df_patient_other_samples$Sample[j]
                
                #path til vcf fil sample j:
                analyseID2<-df_patient_other_samples$Analysis.ID[j]
                sampleID2<-df_patient_other_samples$Sample.ID[j]
                panelID2<-df_patient_other_samples$Panel[j]
                panelID2<-ifelse(panelID2 == "Surveillance v2","Surveillance-v2",panelID2)
                path_sample_mappe2<-paste0(path_AVENIO_results,
                                           "/Plasma-",
                                           analyseID2,
                                           "/",
                                           sampleID2,
                                           "-",
                                           panelID2)
                path_VCF_fil<-paste0(path_sample_mappe2,
                                     "/snv-RocheDefault-",
                                     sampleID2,
                                     ".VCF")
                print(path_VCF_fil)
                #gem VCF fil:
                vcf_fil<-vcfR::read.vcfR(path_VCF_fil)
                
                #hvis varianter fundet, tjekkes disse i bam filen en efter en:
                if (nrow(vcf_fil) == 0) 
                {print(paste0("No variants in ", sample2))}
                else{
                    for (variant in 1:nrow(vcf_fil)){
                        vcf_variant<-vcf_fil[variant]
                        
                        #tjek om variant allerede er fundet i sample1:
                        x<-df_nested$variants[[i]][,1:2]
                        chrom<-vcf_fil[variant]@fix[1]
                        pos<-vcf_fil[variant]@fix[2]
                        
                        if (sum(chrom==x$chrom&pos==x$pos)==0){
                            
                            reads_variant_pos<-AvenioUpdate:::get_reads_at_variant_position_using_VCF_in_deduped(
                                vcf_variant,sampleID,path_sample_mappe)
                            variant_depth_at_pos<-AvenioUpdate:::get_variant_depth_at_variant_position_using_VCF(
                                vcf_variant,reads_variant_pos)
                            liste_variant<-c(c(sample1,sample2),variant_depth_at_pos)
                            #hvis varianten er fundet med mininum 1 read, tilføjes denne:
                            if (liste_variant[4]!="0"){
                                df_varianter_fundet_i_bam<-rbind(df_varianter_fundet_i_bam,
                                                                 liste_variant)
                            }}
                    }}
            }}
    }
    #slet den tomme række
    df_varianter_fundet_i_bam<-df_varianter_fundet_i_bam[-1,]
    return(df_varianter_fundet_i_bam)
}