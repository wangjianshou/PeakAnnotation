#!/PUBLIC/software/RNA/R-3.1.2/R-3.1.2/bin/Rscript
# PeakStartCodon.r gene.gtf *_peaks.narrowPeak...
# argument <- c("gene.gtf", "peak")
argument <- commandArgs(TRUE)
options(stringsAsFactors=FALSE)
options(digits=10)
.libPaths("/ifs/TJPROJ3/HW/wangjianshou/R_lib/3.1")
library(data.table)
library(reshape2)
library(ggplot2)

ori_gtf <- fread(argument[1], header=F, data.table=T, sep="\t")
setnames(ori_gtf, c("chr", "src", "type", "start_pos", "end_pos", "score", "strand", "codon", "info"))
ori_gtf <- ori_gtf[type %in% c("exon", "start_codon")]
ori_gtf[, `:=`(gene_id=gsub(".*gene_id \"(.*?)\";.*", "\\1", info),
               transcript_id=gsub(".*transcript_id \"(.*?)\".*", "\\1", info),
               exon_num=as.integer(gsub(".*exon_number \"(\\d+)\".*", "\\1", info)),
               tr_biotype=gsub(".*transcript_biotype \"(.*?)\".*", "\\1", info),
               gene_biotype=gsub(".*gene_biotype \"(.*?)\".*", "\\1", info),
               info=NULL,
               src=NULL,
               score=NULL,
               codon=NULL
			  )
	   ]
StartCodonTr <- ori_gtf[type=="start_codon"]$transcript_id
ori_gtf <- ori_gtf[transcript_id %in% StartCodonTr]
gene_length <- ori_gtf[type=="exon",
                      .(alength=sum(end_pos-start_pos+1)),
                      by=.(gene_id, transcript_id)
                      ][,
                      .(transcript_id=transcript_id[which.max(alength)], alength=max(alength)),
                      by = .(gene_id)
                      ]
ori_gtf <- ori_gtf[transcript_id %in% gene_length$transcript_id]
ori_gtf[, cls := exon_num[which(type=="start_codon")], by=transcript_id]
gtf <- ori_gtf[type=="exon"][, alength := end_pos-start_pos+1]

gtf[order(exon_num),
    be_length := cumsum(append(alength[-1], c(0, 0), 0)[-(length(alength)+1)]),
    by=transcript_id
     ]
write.table(gtf, file="myExonGTF", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

setkey(gtf, chr, start_pos)
for (chr in unique(gtf$chr)){
    eval(parse(text=paste(chr, " <- gtf[chr==", "\"", chr, "\"", "]", sep="")))
}

peakMock <- fread("Mock_peaks.narrowPeak", select=1:4, header=FALSE, sep="\t", data.table=TRUE)
setnames(peakMock, c("chr", "start_pos", "end_pos", "peak_name"))
peakMock <- peakMock[order(chr, start_pos)][, middle := floor(start_pos+(end_pos-start_pos)/2)]
peakMock <- peakMock[, head(.SD, 1), by=.(chr, middle)]

fun <- function(Vchr, Vpos) {
    lin <- with(get(Vchr), which(Vpos >= start_pos & Vpos < end_pos))
    if(length(lin)==0){
        return(list("NA", 1024L, 10241024L))
        }
    else if(length(lin) == 1){
        return(get(Vchr)[lin,
                        .(transcript_id, cls,
                          if(exon_num==1)
                            ifelse(strand=="+",
                                   as.integer(Vpos-end_pos),
                                   as.integer(start_pos-Vpos)
                                   )
                             else
                              ifelse(strand=="+",
                                     as.integer(Vpos-start_pos+be_length),
                                     as.integer(end_pos-Vpos+be_length)
                                     )
                          )
                        ]
               )
        }
    else if(length(lin) > 1){
        return(get(Vchr)[lin[sample(1:length(lin), size=1)],
                        .(transcript_id, cls,
                          if(exon_num==1)
                            ifelse(strand=="+",
                                   as.integer(Vpos-end_pos),
                                   as.integer(start_pos-Vpos)
                                   )
                          else
                            ifelse(strand=="+",
                                   as.integer(Vpos-start_pos+be_length),
                                   as.integer(end_pos-Vpos+be_length)
                                   )
                          )
                        ]
               )
        }
    else{
        return(list("NA", -1024, -10241024))
        }
}
peakMock[, c("transcript_id", "cls", "distance") := fun(chr, middle), by=.(chr, middle)]
peakMock <- peakMock[transcript_id != "NA"]
write.table(peakMock, file="peakMock", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

peakALKBH3 <- fread("ALKBH3_peaks.narrowPeak", select=1:4, header=FALSE, sep="\t", data.table=TRUE)
setnames(peakALKBH3, c("chr", "start_pos", "end_pos", "peak_name"))
peakALKBH3 <- peakALKBH3[order(chr, start_pos)][, middle := floor(start_pos+(end_pos-start_pos)/2)]
peakALKBH3 <- peakALKBH3[, head(.SD, 1), by=.(chr, middle)]
peakALKBH3[, c("transcript_id", "cls", "distance") := fun(chr, middle), by=.(chr, middle)]
peakALKBH3 <- peakALKBH3[transcript_id != "NA"]
write.table(peakALKBH3, file="peakALKBH3", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

Mock <- peakMock[, .(peak_name, transcript_id, cls, distance)
                ][, clas := ifelse(cls==1, "aaa", ifelse(cls==2, "bbb", "ccc"))
                ][, cls := NULL]
MockAll <- Mock[, .SD][, clas := "ddd"]
Mock <- rbind(Mock, MockAll)



#密度图
MockFig <- ggplot(data=Mock, aes(x=distance, group=clas)) +
           geom_line(aes(y=..count.., color=clas), adjust=0.5, stat="density") +
           labs(title="Mock", x="Distance from 1st splice site(nucleotides)", y="Number of peaks") +
           scale_x_continuous(expand=c(0,0), limits=c(-3000,5000)) +
           scale_color_hue(labels=c(aaa="AUG in 1st exon", bbb="AUG in 2nd exon", ccc="AUG in 3rd+ exon", ddd="All genes")) +
           theme(
           axis.title.x=element_text(hjust=0.5, vjust=0),
           axis.title.y=element_text(hjust=0.6, vjust=0.5),
           legend.title=element_blank(),
           plot.title=element_text(size=15, vjust=1)
           )
ggsave(file="Mock.png", plot=MockFig, type="cairo", dpi=600)

ALKBH3Fig <- ggplot(data=ALKBH3, aes(x=distance, group=clas)) +
           geom_line(aes(y=..count.., color=clas), adjust=0.5, stat="density") +
           labs(title="ALKBH3", x="Distance from 1st splice site(nucleotides)", y="Number of peaks") +
           scale_x_continuous(expand=c(0,0), limits=c(-3000,5000)) +
           scale_color_hue(labels=c(aaa="AUG in 1st exon", bbb="AUG in 2nd exon", ccc="AUG in 3rd+ exon", ddd="All genes")) +
           theme(
           axis.title.x=element_text(hjust=0.5, vjust=0),
           axis.title.y=element_text(hjust=0.6, vjust=0.5),
           legend.title=element_blank(),
           plot.title=element_text(size=15, vjust=1)
           )
ggsave(file="ALKBH3.png", plot=ALKBH3Fig, type="cairo", dpi=600)






























