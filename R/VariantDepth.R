#' @noRd
VariantDepth <- function(Chr, Pos, Base, Bam) {
    base <- substring(as.character(Bam[[1]]$seq),
                      Pos-Bam[[1]]$pos+1,
                      Pos-Bam[[1]]$pos+1)
    SeqVariant<-Bam[[1]]$seq[base==Base]
    depth<-length(SeqVariant)
    depth
}