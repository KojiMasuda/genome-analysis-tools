#ifndef _PARSE_CHR_H_
#define _PARSE_CHR_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>

#include "write_tab.h"

#define LINE_STR_LEN 100000 //char length per line

/*
 * Structure of chr block.
 * This is a link list.
 * For each chr block, binding site(bs) list is linked.
 */
struct chr_block {
  char *chr;
  struct chr_block *next;
  struct bs *bs_list;
  struct sig *sig_list;
  unsigned long bs_nb;
  int bs_init;
  int sig_init;
};

/*
 * Structure of bs.
 * This is a link list.
 * Binding site(bs) list is linked to one chr block.
 */
struct bs {
  unsigned long st;
  unsigned long ed;
  char strand;
  char *line;
  struct bs *next;
  struct bs *prev;
};

/*
 * Structure of sig.
 * This is a link list.
 * sig data list is linked to one chr block.
 */
struct sig {
  unsigned long st;
  unsigned long ed;
  float val;
  struct sig *next;
  struct sig *prev;
};

extern char *ga_header_line;

void ga_parse_chr_bs (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int col_strand, int hf);
void ga_parse_chr_bs_rand (struct chr_block **chr_block_head, struct chr_block *chr_block_head_ori, struct chr_block *chr_table);
void ga_parse_bedgraph (const char *filename, struct chr_block **chr_block_head);
void ga_parse_sepwiggz (const char *filename, struct chr_block **chr_block_head);
void ga_parse_onewiggz (const char *filename, struct chr_block **chr_block_head);
void ga_free_chr_block (struct chr_block **chr_block);
unsigned long ga_count_peaks (struct chr_block *chr_block_head);

#endif
