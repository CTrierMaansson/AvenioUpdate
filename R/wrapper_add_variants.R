#' @noRd
#' @importFrom dplyr distinct `%>%`
wrapper_add_variants<-function(df1,df2){
    #fjern varianter fra df2, der går igen:
    df2<-df2%>%dplyr::distinct(
        sample_1,genomic_position,
        variant_depth,unique_depth,
        .keep_all = TRUE)
    
    #tilføj alle varianter en efter en:
    for (i in 1:nrow(df2)){
        df1<-add_bam_mutation(df = df1,
                              Sample1 = df2$sample_1[i],
                              Sample2 = df2$sample_2[i],
                              GenomicPosition = df2$genomic_position[i],
                              VariantDepth = as.numeric(df2$variant_depth[i]),
                              UniqueDepth = as.numeric(df2$unique_depth[i]))
    }
    return(df1)
}