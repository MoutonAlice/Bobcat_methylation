Here is a tutorial with methylkit followed by genomation:
  
  https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#4_annotating_differentially_methylated_bases_or_regions

Below my code: 
  
  # to transform the gtf file into the genomation format we use two tools available in UCSC (but you can also download the right format directly from UCSC, te cat genome is not updated so I had to find another solution)
  
  ./gtfToGenePred Feliscatus.gtf FeliscatusGenepred
./genePredToBed FeliscatusGenepred Feliscatus12.bed

gene.obj=readTranscriptFeatures("/u/flashscratch/flashscratch1/a/amouton/genomes/cat/Feliscatus12.bed") # path to your bed file (computer or hoffman)
macau<-read.table("macau_genomation.txt",header=T)
#   chr     start       end   beta        qval
# 1  A1   5383180   5383180 -20.70 0.016433073
# 2  A1  66253437  66253437 -15.18 0.067072388
# 3  A1 112874794 112874794 -21.01 0.000086600
# 4  A1 112899729 112899729  19.96 0.050897578
# 5  A1 116308397 116308397 -20.40 0.005621199
# 6  A1 226485133 226485133  18.32 0.028136964

### annotate differentially methylated CpGs with promoter/exon/intron using annotation data
annotateWithGeneParts(as(macau,"GRanges"),gene.obj)

### getAssociationWithTSS function from genomation package.
#After getting the annotation of differentially methylated regions, we can get the # distance to TSS and nearest gene name using the  : 
diffAnn=annotateWithGeneParts(as(macau,"GRanges"),gene.obj)
annotation<-getAssociationWithTSS(diffAnn)

# membership of each site (promoter intro exon)
member<-getMembers(diffAnn)
write.table(annotation,file="annotation_macau.txt",sep="\t")
write.table(member,file="member_macau.txt",sep="\t")

# combine in one file everything
merged.table<-cbind(macau[annotation$target.row],annotation)

#### plot methylation pattern
png("plotmethylation_genomation_macau.png")
plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")
dev.off()

## read the CpG island annotation and annotate our differentially methylated bases/regions with them.

# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank("/u/flashscratch/flashscratch1/a/amouton/genomes/cat/Feliscatus12.bed",feature.flank.name=c("CpGi","shores","enhancers"))#

#### convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(macau,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
member<-getMembers(diffCpGann)
annotation_cpg<-getAssociationWithTSS(diffCpGann)
write.table(annotation,file="annotation_macau.txt",sep="\t")
write.table(member,file="member_cpgisland_macau.txt",sep="\t")


png("plotmethylation_cpgisland.png")
plotTargetAnnotation(diffCpGann,precedence=TRUE,
                     main="differential methylation annotation cpgisland")
dev.off()


## Match gene name and transcript to make sense
gtf1 = rtracklayer::import("./Singlesite_filtering5/18samples/genomation/Feliscatus.gtf") # path to your gtf file (computer or hoffman)
gtf1 = as.data.frame(gtf1[gtf1$type == "gene",])
write.table(gtf1[,c("gene_id", "gene_name", "gene_biotype","strand","transript_id")], file = "gene_names.txt", col.names = T, row.names = F, sep = "\t", quote = F)

Gene_name <-read.table("gene_names.txt",header=T)
head(Gene_name)
#order the  id
Gene_name<-Gene_name[order(Gene_name$transcript_id),]
match<-merge(Gene_name,annotation,by.x="transcript_id",by.y="feature.name")
dim(match) # 
match0<-match[!duplicated(match$gene_name),]
dim(match0) # 
write.table(match0,file="gene_name_methylationsite.txt",sep="\t")
