genome-analysis-tools
========
These tools are created for analyzing the peaks or read distributions (bedgraph, wiggle format) derived from next-generation sequencing such as ChIP-seq or RNA-seq.

The included tools are:
* `ga_overlap`: check the overlapping and return the overlapping, non-overlapping, and original file with ov/nonov flags.
* `ga_reads_summit`: calculate read distributions around summit or specific position such as TSS (Transcription Start Site).
* `ga_reads_gene`: calculate read distributions around genes.
* `ga_intersummit`: calculate inter-summit distances.
* `ga_reads_region`: calculate read amounts in the specific regions such as peak, up-stream regions.
* `ga_flanking`: pick up the regions (genes) which flank peaks/summits.
* `ga_nuc_summit`: calculte nucleotide content distributions around summits. User-specified motif distribution can be calculated.

Requirements
========
* [GCC compiler](http://gcc.gnu.org/)
