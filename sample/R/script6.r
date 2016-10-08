source("sample/R/heatmap_f.r") ##reading tiny R script
cols <- colorRampPalette(c("white","darkred")) ##setting colors from white to darkred

df <- read.delim ("sample/chip1/chip1_divided_input1_around_peak1_halfwid1000winsize300step10_all.txt") ##reading output file of ga_reads_summit_all
df.sort <- df[order(apply(df[,91:111], MARGIN=1, mean)),] ##sorting by means of win for -100 bp to +100 bp (win for center, 0 bp, is element of 101)
df.sort.norm.all <- return_mat_norm_all(df.sort, 0.99)

png("sample/chip1/chip1_divided_input1_around_peak1_halfwid1000winsize300step10_all.png",height=675)
image(seq(from=-1000, to=1000, by = 10),1:nrow(df.sort.norm.all),as.matrix(t(df.sort.norm.all)), col= cols(100), main="Enrichment of chip1 around peak1\nsorted by -100 to 100 bp enrichment", xlab="Distance from peak1 summits (bp)", ylab="peak1", xlim=c(-1000, 1300), bty="n") ; abline(v=0)
rect(1100, seq(from=10, to=109, by=1), 1200, seq(from=11, to=110,by=1),col=cols(100),border="transparent")
text(x=1250, y=10, "0")
text(x=1250, y=109, "1")
dev.off()
