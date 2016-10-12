df <- read.delim ("sample/peak1_nuc_hw1000_win100_step10_AT.txt") ##reading output file

png("sample/peak1_nuc_hw1000_win100_step10_AT.png")
plot (seq(from = -1000, to = 1000, by = 10), apply(df, MARGIN=2, mean) * 100, col="red", type="l", xlab="Distance from peak1 summit (bp)", ylab="AT content (%)", main="AT content around peak1", lwd=2)
abline(v=0, lty=2)
dev.off()
