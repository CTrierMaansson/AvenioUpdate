#' @noRd
add_variants_from_bam_files<- function(df){
    #tjek bamtiler for snv og indels
    snv_fundet_i_bam<-wrapper_tjek_bam_snv(df)
    indel_findet_i_bam<-wrapper_tjek_bam_indel(df)
    
    #tilføj varianter til df:
    if(nrow(snv_fundet_i_bam)>0){
        df<-wrapper_add_variants(df,snv_fundet_i_bam)}
    
    if(nrow(indel_findet_i_bam)>0){
        df<-wrapper_add_variants(df,indel_findet_i_bam)}
    
    return(df)
}