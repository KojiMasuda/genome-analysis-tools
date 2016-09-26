genome-analysis-tools
========
These tools are created for analyzing the peaks or read distributions (bedgraph, wiggle format) derived from next-generation sequencing such as ChIP-seq or RNA-seq.

The included tools are:
* `ga_overlap`: checks the overlapping and return the overlapping, non-overlapping, and original file with ov/nonov flags.
* `ga_reads_summit`: calculates the average read distribution around summit or specific position such as TSS (Transcription Start Site).
* `ga_reads_summit_all`: calculates read distributions around ALL summits or specific positions such as TSS (Transcription Start Site).
* `ga_calc_dist`: calculates the distance between two peak sets or inter-summit distance.
* `ga_reads_region`: calculates read amounts in the specific regions such as peak, up-stream regions.
* `ga_nuc_region`: calcultes nucleotide content in the specified regions.
* `ga_nuc_summit`: calcultes nucleotide content distributions around summits.
* `ga_deltaG`: makes the wiggle file of the free energy difference between the duplex and single-strand states from fasta file.
* `ga_RPKM`: calculates the expression levels as RPKM.
* `ga_reads_gene`: calculates read distributions around genes (not supported yet...).
* `ga_flanking`: picks up the regions (genes) which flank peaks/summits (not supported yet...).

Requirements
========
* [GCC compiler](http://gcc.gnu.org/)
* [GNU C Library(glibc)](http://www.gnu.org/software/libc/)
* [zlib](http://www.zlib.net/)
