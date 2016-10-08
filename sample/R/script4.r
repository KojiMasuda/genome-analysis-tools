df <- read.delim ("sample/nuc.wig_around_anno_halfwid1000winsize25step10.txt")

png("sample/nuc.wig_around_anno_halfwid1000winsize25step10.png")
plot (df$relative_pos, df$smt_mean, col="black", ylim=c(-1,0.35), main="Nucleosome distribution around TSS", ylab="Nucleosome occupancy (z-score)", xlab="Distance from TSS (bp)", type="l", lwd=2)
lines (df$relative_pos, df$CI95.00percent_U, col="black", lty=2)
lines (df$relative_pos, df$CI95.00percent_L, col="black", lty=2)
abline (v = 0, h = 0, lty=2)
dev.off()
