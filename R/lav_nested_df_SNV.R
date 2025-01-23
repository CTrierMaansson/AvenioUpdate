#' @noRd
#' @importFrom dplyr tibble `%>%` distinct filter nest_join
#' @importFrom stringr str_split_i str_sub
lav_nested_df_SNV<-function(df){
    df1<-dplyr::tibble(
        Sample=df$sample_index,
        Patient=stringr::str_split_i(df$sample_index,"_",2),
        Sample.ID=df$Sample.ID,
        Subject.ID=df$Subject.ID,
        Date.of.Sample.Collection=df$Date.of.Sample.Collection,
        Analysis.Name=df$Analysis.Name,
        Analysis.ID=df$Analysis.ID,
        Panel=df$Panel)
    df1<-dplyr::distinct(df1)
    
    #kun SNV'er:
    df<-df%>%dplyr::filter(Mutation.Class=="SNV")
    
    df2<-dplyr::tibble(
        Sample=df$sample_index,
        chrom=df$Genomic.Position%>%stringr::str_split_i(pattern=":",i=1),
        pos=df$Genomic.Position%>%stringr::str_split_i(pattern=":",i=2),
        gene=df$Gene,
        ref=df$Coding.Change%>%
            stringr::str_split_i(pattern=">",i=1)%>%
            stringr::str_sub(-1),
        alt=df$Coding.Change%>%
            stringr::str_split_i(pattern=">",i=2)%>%
            stringr::str_sub(1,1))
    df2[is.na(df2)]<-"NA"
    df3<-dplyr::nest_join(df1,df2,name="variants")
    return(df3)
}
