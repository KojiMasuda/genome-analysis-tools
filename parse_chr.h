#ifndef _PARSE_CHR_H_
#define _PARSE_CHR_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

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
};

/*
 * Structure of bs.
 * This is a link list.
 * Binding site(bs) list is linked to one chr block.
 */
struct bs {
  int st;
  int ed;
  char *line;
//  char *ov;
  struct bs *next;
//  struct bs *prev;
};

extern char *ga_header_line;

void ga_parse_chr_bs (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int hf);
void ga_free_chr_block (struct chr_block **chr_block);

#endif
