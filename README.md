genome-analysis-tools
========
These tools are created for analyzing the peaks or read distributions (bedgraph, wiggle format) derived from next-generation sequencing such as ChIP-seq or RNA-seq.

The included tools are:
* `ga_overlap`: check the overlapping and return the overlapping, non-overlapping, and original file with ov/nonov flags.
* `ga_reads_summit`: calculate the average read distribution around summit or specific position such as TSS (Transcription Start Site).
* `ga_reads_summit_all`: calculate read distributions around ALL summits or specific positions such as TSS (Transcription Start Site).
* `ga_calc_dist`: calculate the distance between two peak sets or inter-summit distance.
* `ga_reads_gene`: calculate read distributions around genes (not supported yet...).
* `ga_reads_region`: calculate read amounts in the specific regions such as peak, up-stream regions (not supported yet...).
* `ga_flanking`: pick up the regions (genes) which flank peaks/summits (not supported yet...).
* `ga_nuc_summit`: calculte nucleotide content distributions around summits. User-specified motif distribution can be calculated (not supported yet...).

Requirements
========
* [GCC compiler](http://gcc.gnu.org/)
* [GNU C Library(glibc)](http://www.gnu.org/software/libc/)
* [zlib](http://www.zlib.net/)
