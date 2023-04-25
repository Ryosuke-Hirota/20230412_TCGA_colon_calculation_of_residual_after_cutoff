# This script is to calculate miRNA biogenesis efficiency and draw plot about correlation between expression levels of miRNA and transcirpt
# miRNA biogenesis efficiency is calculated after cutoff for log2(miRNA level)ÅÜ 0 and sample without expression ÅÖ 220 
# 2023/4/12 made

# activate packages
library(stringr)
library(ggplot2)

# function for calculating outliers
cal.outlier <-function(x){
  q <-as.numeric(quantile(x))
  iqr <-IQR(x)
  outlier1 <-q[2]-iqr*1.5
  outlier2 <-q[4]+iqr*1.5
  outliers <-append(outlier1,outlier2)
  return(outliers)
}

# import correspondence table between TCGA colon transcriptome bam and TCGA colon miRNA quantification file
# this table is located at "https://github.com/Ryosuke-Hirota/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data"
setwd("C:/Rdata/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data")
cor.table <-read.table("correspondence_table_between_TCGA_colon_transcriptome_bam_and_miRNA_qunatification_file.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of transcripts that intersect with miRNAs in gencode v36
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata")
primir.list <-read.table("TCGA_hg38_transcript_intersect_with_miRNA.txt",sep="\t",header = F,stringsAsFactors = F)
primir.list[,2] <-primir.list[,2]-1
primir.list <-primir.list[primir.list[,2]!=primir.list[,8]&primir.list[,3]!=primir.list[,9],]
primir.list <-primir.list[,c(4,11,10)]
transcripts <-unique(primir.list[,3])

# import data about mean of transcript level or miRNA level
# this data is located at "https://github.com/Ryosuke-Hirota/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_l"
setwd("C:/Rdata/20230215_TCGA_colon_percentage_of_combinations_with_sig_posi_correlation_upon_cutoff_of_expression_level")
mean.df <-read.table("TCGA_colon_table_about_number_of_sample_without_expression.txt",sep="\t",header=T,stringsAsFactors = F)

# cutoff (log2[miRNA exp]ÅÜ0, number of sample without expression ÅÖ 220)
# 190 combinations with significant positive correlations
mean.df <-mean.df[mean.df[,6]>=0,]
mean.df <-mean.df[mean.df[,8]<=220,]
nrow(mean.df[mean.df[,3]>0&mean.df[,4]<0.05,]) 

# list TCGA colon transcript quantification files
# these files are located at "\\fsw-q02\okamura-lab\20221006_TCGA_colon_salmon_quant_transcriptome"
setwd("C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification")
transcript.quant <-list.files(path = "C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification",pattern = ".txt")
t.file.id <-gsub("_quant.txt","",transcript.quant)

# make table summarized TCGA colon transcript quantification
for (i in 1:length(transcript.quant)){
  transcript.file <-read.table(transcript.quant[i],sep = "\t",header = T,stringsAsFactors = F)
  t <-match(transcripts,transcript.file[,1])
  transcript.file <-transcript.file[t,]
  transcript.file <-transcript.file[,c(1,4)]
  colnames(transcript.file)[2] <-transcript.quant[i]
  if(i==1){
    transcript.quant.table <-transcript.file
  }else{
    transcript.quant.table <-merge(transcript.quant.table,transcript.file,by="Name")
  }}

# import miRNA quantification table
# this table made from \\fsw-q02\okamura-lab\20221006_TCGA_colon_salmon_quant_transcriptome""
setwd("C:/Rdata/20230105_TCGA_colon_miRNA_quantification")
miRNA.quant.table <-read.table("table_of_TCGA_colon_miRNA_quantifications.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# make new directory
setwd("C:/Rdata")
dir.create("20230412_TCGA_colon_calculation_of_residual_after_cutoff")
setwd("C:/Rdata/20230412_TCGA_colon_calculation_of_residual_after_cutoff")
dir.create("plot")
dir.create("residual")

# investigate correlation and residual between expression level of transcript and miRNA, draw plot and make summary
for (i in 1:nrow(mean.df)){
  # extract expression level of a certain miRNA and edit dataframe 
  miRNA.df <-miRNA.quant.table[miRNA.quant.table[,1]==mean.df[i,1],]
  miRNA.df <-as.data.frame(t(miRNA.df),stringsAsFactors = F)
  m.cor <-match(rownames(miRNA.df),cor.table[,5])
  miRNA.df[,2] <-rownames(miRNA.df)
  miRNA.df[,3] <-cor.table[m.cor,2]
  colnames(miRNA.df) <-c(mean.df[i,1],"miRNA_file_name","transcript_file_name")
  rownames(miRNA.df) <-NULL
  miRNA.df <-miRNA.df[-1,]
  
  # extract expression level of a certain transcript and edit dataframe
  transcript <-primir.list[primir.list[,1]==mean.df[i,1]&primir.list[,2]==mean.df[i,2],3]
  t <-match(transcript,transcript.quant.table[,1])
  transcript.df <-transcript.quant.table[t,]
  transcript.df <-transcript.df[,-1]
  transcript.df[1,] <-apply(transcript.df, 2, sum)
  transcript.df <-transcript.df[c(-2,-3),] 
  transcript.df <-as.data.frame(t(transcript.df),stringsAsFactors = F)
  t.cor <-match(t.file.id,cor.table[,1])
  transcript.df[,2] <-cor.table[t.cor,2]
  colnames(transcript.df) <-c(mean.df[i,2],"transcript_file_name")
  rownames(transcript.df) <-NULL
  
  # merge dataframe about expression levels of miRNA and transcript
  mt.df <-merge(miRNA.df,transcript.df,by="transcript_file_name")
  mt.df <-subset(mt.df,!is.na(mt.df[,1]))
  mt.df <-mt.df[,c(4,2,1,3)]
  mt.df[,2] <-as.numeric(mt.df[,2])
  
  # remove samples without expression levels of transcript and miRNA
  ex.mt.df <-mt.df[mt.df[,1]!=0&mt.df[,2]!=0,]
  
  # normalize with log2 (after excluding samples without expression) 
  ex.mt.df[,1] <-log2(ex.mt.df[,1])
  ex.mt.df[,2] <-log2(ex.mt.df[,2])
  
  # calculate outliers (outlier is Q1-1.5*IQR or Q3+1.5*IQR)
  t.outlier <-cal.outlier(ex.mt.df[,1])
  m.outlier <-cal.outlier(ex.mt.df[,2])
  
  # remove outliers
  ex.mt.df <-ex.mt.df[ex.mt.df[,1]>t.outlier[1]&ex.mt.df[,1]<t.outlier[2],]
  ex.mt.df <-ex.mt.df[ex.mt.df[,2]>m.outlier[1]&ex.mt.df[,2]<m.outlier[2],]
  
  # make empty summary of residual
  residual.sm <-as.data.frame(matrix(nrow = nrow(ex.mt.df),ncol = 3))
  colnames(residual.sm) <-c("residual","transcriptome_bam","miRNA_quantification")
  
  # calculate correlation coefficient
  r <-try(cor.test(ex.mt.df[,1],ex.mt.df[,2],method="pearson"),silent = T)
  
  if(class(r)!="try-error"){
    # perform linear regression analysis and calculate residual
    l <-lm(ex.mt.df[,2]~ex.mt.df[,1],data = ex.mt.df)
    res <-residuals(l)
    
    # draw plot about correlation between expression levels of transcript and miRNA
    p <-ggplot(data=ex.mt.df,aes(x=ex.mt.df[,1],y=ex.mt.df[,2]))+
      geom_point(color="blue")+
      geom_smooth(data=ex.mt.df,mapping = aes(x=ex.mt.df[,1],y=ex.mt.df[,2]),method="lm",formula='y~x',se=FALSE,colour="black",linewidth=0.5)+
      labs(title=paste0("R =",signif(r$estimate,3),", p = ",signif(r$p.value,3),", n = ",nrow(ex.mt.df)),x=paste0("log2(",mean.df[i,2],")"),
           y=paste0("log2(",mean.df[i,1],")"))+ 
      theme_bw()+
      theme(legend.background = element_rect(fill = "white", colour = "black"))
    
    # save plot about correlation between expression levels of transcript and miRNA
    setwd("C:/Rdata/20230412_TCGA_colon_calculation_of_residual_after_cutoff/plot")
    ggsave(filename=paste0("plot_of_correlation_between_",mean.df[i,1],"_and_",mean.df[i,2],".pdf"),plot = p)
    
    # write summary of residual
    residual.sm[,1] <-res
    residual.sm[,2] <-ex.mt.df[,3]
    residual.sm[,3] <-ex.mt.df[,4]
  
    # output summary of residual
    setwd("C:/Rdata/20230412_TCGA_colon_calculation_of_residual_after_cutoff/residual")
    write.table(residual.sm,paste0("table_of_residual_about_",mean.df[i,1],"_vs_",mean.df[i,2],".txt"),sep="\t",row.names = F,quote = F)
}}

