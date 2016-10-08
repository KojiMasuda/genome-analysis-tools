df <- read.delim ("sample/chip1/chip1_around_peak1_halfwid1000winsize25step10.txt")

png("sample/chip1/chip1_around_peak1_halfwid1000winsize25step10.png")
plot (df$relative_pos, df$smt_mean, col="blue", ylim=c(4,15), main="Read distribution of chip1 around peak1", ylab="Read density", xlab="Distance from peak1 (bp)", type="l", lwd=2)
lines (df$relative_pos, df$CI95.00percent_U, col="blue", lty=2)
lines (df$relative_pos, df$CI95.00percent_L, col="blue", lty=2)
dev.off()
