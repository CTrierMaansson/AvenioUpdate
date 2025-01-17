#' @noRd
#' @importFrom dplyr tibble
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr nest_join
#' @importFrom stringr str_split_i
#' @importFrom stringr str_sub
#' @importFrom stringr str_glue
#' @importFrom string str_extract
lav_nested_df_indel<-function(df){
    df1<-dplyr::tibble(
        Sample=df$sample_index,
        Patient=str_split_i(df$sample_index,"_",2),
        Sample.ID=df$Sample.ID,
        Subject.ID=df$Subject.ID,
        Date.of.Sample.Collection=df$Date.of.Sample.Collection,
        Analysis.Name=df$Analysis.Name,
        Analysis.ID=df$Analysis.ID)
    df1<-dplyr::distinct(df1)
    #kun INDEL:
    df<-df%>%dplyr::filter(Mutation.Class=="INDEL")
    #lav liste med type:
    type_list<-sapply(df$Coding.Change, function(x){
        if (!is.na(stringr::str_extract(x,"del"))){
            antal<-as.character(nchar(stringr::str_split_i(x,"del",i=-1)))
            return(stringr::str_glue("{antal}D"))
        }else if (!is.na(stringr::str_extract(x,"dup"))){
            antal<-as.character(nchar(stringr::str_split_i(x,"dup",i=-1)))
            return(stringr::str_glue("{antal}I"))
        }
    })
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
            stringr::str_sub(1,1),
        type=type_list
    )
    df2[is.na(df2)]<-"NA"
    df3<-dplyr::nest_join(df1,df2,name="variants")
    return(df3)
}