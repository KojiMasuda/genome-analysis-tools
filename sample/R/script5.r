source("sample/R/heatmap_f.r") ##reading tiny R script
cols <- colorRampPalette(c("white","darkblue")) ##setting colors from white to darkblue

df <- read.delim ("sample/chip2/chip2_around_anno_mm10_halfwid5000winsize500step100_all.txt") ##reading output file of ga_reads_summit_all
set.seed(127); df.r <- df[sample(1:nrow(df), 2500),] ##randomly picking up 2500 genes
ot <- read.delim ("sample/order_table.txt", header=F) ##order table for sorting...

df.r.sort <- df.r[ot[,1],] ##sorting

df.r.sort.norm <- return_mat_norm_all(df.r.sort, 0.99) ##normalization considering all peaks, where top 99% and more was set to 1
png("sample/chip2/chip2_around_anno_mm10_halfwid5000winsize500step100_all.png",height=1000)
image(seq(from=-5000, to=5000, by = 100),1:nrow(df.r.sort.norm),as.matrix(t(df.r.sort.norm)), col= cols(100), main="Read density of chip2 around TSS", xlab="Distance from TSS (bp)", ylab="Genes", xlim=c(-5000,5900), bty="n") ; abline(v=0)
rect(5200, seq(from=10, to=1000, by=10), 5700, seq(from=20, to=1010,by=10),col=cols(100),border="transparent")
text(x=5800, y=10, "0")
text(x=5800, y=1000, "1")
dev.off()
