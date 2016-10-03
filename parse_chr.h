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
  struct chr_block *tail;
  struct bs *bs_list;
  struct sig *sig_list;
  struct ref *ref_list;
  unsigned long bs_nb;
  int bs_init;
  int sig_init;
  int ref_init;
};

/*
 * Structure of chr block fa.
 * This is a link list.
 * For each chr block, base letters are stored.
 */
struct chr_block_fa {
  char *chr;
  char *letter;
  struct chr_block_fa *next;
  struct chr_block_fa *tail;
  unsigned long letter_len;
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

/*
 * Structure of ref.
 * This is a link list.
 * reference list is linked to one chr block.
 * *ex_st, *ex_ed, *line must be freed.
 */
struct ref {
  unsigned long st;
  unsigned long ed;
  char strand;
  char *ex_st;
  char *ex_ed;
  char *line;
  char *gene;
  char *rm_ex_st;
  char *rm_ex_ed;
  char *ov_gene;
  double ov_prop;
  struct ref *next;
  struct ref *tail;
  struct ref *prev;
};

extern char *ga_header_line;

void ga_parse_chr_bs (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int col_strand, int hf);
void ga_parse_chr_bs_rand (struct chr_block **chr_block_head, struct chr_block *chr_block_head_ori, struct chr_block *chr_table, int hw);
int ga_parse_chr_ref (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int col_strand, int col_ex_st, int col_ex_ed, int col_gene, int hf);
int ga_parse_chr_fa (const char *filename, struct chr_block_fa **chr_block_head, struct chr_block *chr_block_head_gt);
void ga_parse_bedgraph (const char *filename, struct chr_block **chr_block_head);
void ga_parse_sepwiggz (const char *filename, struct chr_block **chr_block_head);
void ga_parse_onewiggz (const char *filename, struct chr_block **chr_block_head);
void ga_free_chr_block (struct chr_block **chr_block);
void ga_free_chr_block_fa (struct chr_block_fa **chr_block);
unsigned long ga_count_peaks (struct chr_block *chr_block_head);

#endif
