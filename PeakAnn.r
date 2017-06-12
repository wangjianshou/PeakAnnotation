#!/PUBLIC/software/RNA/R-3.1.2/R-3.1.2/bin/Rscript
# PeakAnno.R gene.gtf *_peaks.narrowPeak
# argument <- c("gene.gtf", "peak.head")
options(stringsAsFactors=FALSE)
options(digits=10)
.libPaths("/ifs/TJPROJ3/HW/wangjianshou/R_lib/3.1")
library(data.table)
library(reshape2)
library(ggplot2)
argument <- commandArgs(TRUE)
ori_gtf <- fread(argument[1], header=F, data.table=T,sep="\t")
ori_gtf <- ori_gtf[V3 %in% c("exon", "CDS", "five_prime_utr", "three_prime_utr")]
ori_gtf[, `:=`(V1=paste("chr", V1, sep=""),
               gene_id=gsub(".*gene_id \"(\\w+)\";.*", "\\1", V9),
               transcript_id=gsub(".*transcript_id \"(\\w+)\".*", "\\1", V9),
               exon_num=as.integer(gsub(".*exon_number \"(\\d+)\".*", "\\1", V9)),
               V9=NULL,
               V2=NULL,
               V6=NULL
			   )
	   ]
gene_length <- ori_gtf[V3=="exon",
                      .(alength=sum(V5-V4+1)),
                      by=.(gene_id, transcript_id)
                      ][,
                      .(transcript_id=transcript_id[which.max(alength)], alength=max(alength)),
                      by = .(gene_id)
                      ]
gtf <- ori_gtf[transcript_id %in% gene_length$transcript_id
              ][
              V3 %in% c("CDS", "five_prime_utr", "three_prime_utr")
              ][
              V3 %in% c("five_prime_utr", "three_prime_utr"),
              exon_num := 0
              ]
setkey(gtf, V1, V4)
gtf[V7=="+" & exon_num==0,
   exon_num := 1:.N,
   by=.(V3, transcript_id)
   ][V7=="-" & exon_num==0,
    exon_num := .N:1,
    by=.(V3, transcript_id)
    ]
gtf[,
    alength := V5-V4+1,
    ][,
     sum_length := sum(alength),
     by = .(V3, transcript_id)
     ]

gtf[order(exon_num),
   be_length := cumsum(append(alength[order(exon_num)], 0, 0)[-(length(alength)+1)]),
   by = .(V3, transcript_id)
   ]
write.table(gtf, file="myGTF", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


for (chr in unique(gtf$V1)){
    eval(parse(text=paste(chr, " <- gtf[V1==", "\"", chr, "\"", "]", sep="")))
}

fun <- function(Vx, x) {
    lin <- with(get(Vx), which(x >= V4 & x < V5))
    if(length(lin)==0){
        return(list("chrNA", -1))
        }else{
        return(get(Vx)[lin,
               .(V3, ifelse(V7=="+", (x-V4+be_length)/sum_length,(V5-x+be_length)/sum_length))
               ])
    }
}


peak1 <- fread(argument[2],select=1:4, header=FALSE, sep="\t", data.table=TRUE)
peak1 <- peak1[order(V1,V2)][, middle := floor(V2+(V3-V2)/2)]
peak1 <- peak1[, head(.SD, 1), by=.(V1, middle)]
peak1[, c("region", "percent") := fun(V1, middle), by=.(V1, middle)]
peak2 <- fread(argument[3],select=1:4, header=FALSE, sep="\t", data.table=TRUE)
peak2 <- peak2[order(V1,V2)][, middle := floor(V2+(V3-V2)/2)]
peak2 <- peak2[, head(.SD, 1), by=.(V1, middle)]
peak2[, c("region", "percent") := fun(V1, middle), by=.(V1, middle)]


if(FALSE)
{
    write.table(peak,
            file="Mock_peak_annotation",
            quote=FALSE,
            sep="\t",
            row.names=FALSE,
            col.name=FALSE
            )
    write.table(peak[percent > 0],
            file="Mock_peak_annotation_noNA",
            quote=FALSE,
            sep="\t",
            row.names=FALSE,
            col.name=FALSE
            )
}

peak1_ann_noNA <- peak1[percent > 0]
peak1_pie <- peak1_ann_noNA[, .(num_peak=.N), by=region]
peak1_pie[, prop:=round(num_peak/nrow(peak1_ann_noNA)*100, 2), by=region]
peak1_pie[, labs := switch(region,
                          three_prime_utr=paste("3'UTR\n", prop, "%", sep=""),
                          five_prime_utr=paste("5'UTR\n", prop, "%", sep=""),
                          CDS=paste("CDS\n", prop, "%",sep="")),
        by=region]
#write.table(peak1_pie, file="peak_pie_table", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
peak2_ann_noNA <- peak2[percent > 0]
peak2_pie <- peak2_ann_noNA[, .(num_peak=.N), by=region]
peak2_pie[, prop:=round(num_peak/nrow(peak2_ann_noNA)*100, 2), by=region]
peak2_pie[, labs := switch(region,
                          three_prime_utr=paste("3'UTR\n", prop, "%", sep=""),
                          five_prime_utr=paste("5'UTR\n", prop, "%", sep=""),
                          CDS=paste("CDS\n", prop, "%",sep="")),
        by=region]
#write.table(peak2_pie, file="peak_pie_table", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#饼图
pdf("pie.pdf", width=17.1, height=10)
split.screen(c(1,2))
screen(1)
par(mar=c(5, 2, 2, 1))
pie(peak1_pie$num_peak,
    labels=peak1_pie$labs,
    col=hsv(h=c(6.5,7,8)/12),
    border="white",
    radius=0.5,
    cex=1.2
    )
mtext("Peak1", adj=0.5, side=1, cex=2,line=-5)
screen(2)
par(mar=c(5, 2, 2, 1))
pie(peak2_pie$num_peak,
    labels=peak2_pie$labs,
    col=hsv(h=c(6.5,7,8)/12),
    border="white",
    radius=0.5,
    cex=1.2
    )
mtext("Peak2", adj=0.5, side=1, cex=2,line=-5)
dev.off()

png("pie.png", width=17.1, height=10, units="cm", res=300, type="cairo")
split.screen(c(1,2))
screen(1)
par(mar=c(4, 2, 2, 1))
pie(peak1_pie$num_peak,
    labels=peak1_pie$labs,
    col=hsv(h=c(6.5,7,8)/12),
    border="white",
    radius=0.5,
    cex=0.7
    )
mtext("Peak1", adj=0.5, side=1, cex=1,line=-1)
screen(2)
par(mar=c(4, 2, 2, 1))
pie(peak1_pie$num_peak,
    labels=peak1_pie$labs,
    col=hsv(h=c(6.5,7,8)/12),
    border="white",
    radius=0.5,
    cex=0.7
    )
mtext("Peak2", adj=0.5, side=1, cex=1,line=-1)
dev.off()


standard_len <- gtf[, .(region_len=sum(alength)),
    by=.(V3, transcript_id)
    ][, .(aver_len=sum(region_len)/(.N)),
    by=V3
    ]
standard_len[, len_factor := aver_len/standard_len[V3=="five_prime_utr", aver_len]]
len_fact <- standard_len$len_factor
names(len_fact) <- standard_len$V3


peak1[region=="five_prime_utr", pos_stand := percent]
peak1[region=="CDS",
      pos_stand := percent*len_fact["CDS"] + len_fact["five_prime_utr"]
     ]
peak1[region=="three_prime_utr",
      pos_stand := percent*len_fact["three_prime_utr"] + 
                   len_fact["five_prime_utr"] +
                   len_fact["CDS"]
     ]
peak2[region=="five_prime_utr", pos_stand := percent]
peak2[region=="CDS",
      pos_stand := percent*len_fact["CDS"] + len_fact["five_prime_utr"]
     ]
peak2[region=="three_prime_utr",
      pos_stand := percent*len_fact["three_prime_utr"] + 
                   len_fact["five_prime_utr"] +
                   len_fact["CDS"]
     ]

vline <- cumsum(len_fact[c("five_prime_utr", "CDS", "three_prime_utr")])
peak1[, sample := "Sample1"]
peak2[, sample := "Sample2"]
figdata <- rbind(peak1[, .(sample, pos_stand)], peak2[, .(sample, pos_stand)])


#密度图
peak_distribute_fig <- ggplot(data=figdata, aes(x=pos_stand, group=sample)) +
                              geom_line(aes(color=sample), stat="density") +
                              scale_color_manual(values=hsv(h=c(1,7)/12)) +
                              labs(x="Position in normalized transcript", y="Percent of peaks") +
                              scale_x_continuous(expand=c(0, 0)) +
                              ylim(-0.010, 0.18) +
                              theme(axis.line.x=element_blank(),
                                    panel.grid=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.title.x=element_text(size=15),
                                    axis.title.y=element_text(hjust=0.75, vjust=0.3, size=15)
                                    ) +
                              geom_rect(aes(xmin=0,xmax=vline[1],ymin=-0.008,ymax=-0.004),fill="green") +
                              geom_rect(aes(xmin=vline[1], xmax=vline[2], ymin=-0.010, ymax=-0.002), fill="blue") +
                              geom_rect(aes(xmin=vline[2], xmax=vline[3], ymin=-0.008, ymax=-0.004), fill="yellow") +
                              geom_text(aes(x=vline[1]/2, y=-0.010, label="5'UTR"), vjust=2, size=3) +
                              geom_text(aes(x=(vline[1]+vline[2])/2, y=-0.010, label="CDS"), vjust=2, size=3) +
                              geom_text(aes(x=(vline[2]+vline[3])/2, y=-0.010, label="3'UTR"), vjust=1.5, size=3) +
                              geom_linerange(aes(x=vline[1], ymin=0, ymax=0.17), size=0.5, linetype="dotted") +
                              geom_linerange(aes(x=vline[2], ymin=0, ymax=0.17), size=0.5, linetype="dotted") +
                              geom_text(aes(x=vline[1], y=0.175, label="Start")) +
                              geom_text(aes(x=vline[2], y=0.175, label="Stop"))
ggsave(file="Peak_Distribution.png", plot=peak_distribute_fig, type="cairo", dpi=300)
ggsave(file="Peak_Distribution.samll.pdf", plot=peak_distribute_fig, width=10, height=10, units="cm")



