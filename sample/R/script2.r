df <- read.delim ("sample/chip1/chip1_divided_input1_around_peak1_halfwid1000winsize25step10.txt")

png("sample/chip1/chip1_divided_input1_around_peak1_halfwid1000winsize25step10.png")
plot (df$relative_pos, df$smt_mean, col="blue", ylim=c(1,13), main="Enrichment of chip1 around peak1", ylab="Enrichment of chip1", xlab="Distance from peak1 (bp)", type="l", lwd=2)
lines (df$relative_pos, df$CI95.00percent_U, col="blue", lty=2)
lines (df$relative_pos, df$CI95.00percent_L, col="blue", lty=2)
dev.off()
