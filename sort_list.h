#ifndef _SORT_LIST_H_
#define _SORT_LIST_H_

#include "parse_chr.h"

#include <string.h>
#include <stdio.h>

struct chr_block *ga_mergesort_chr(struct chr_block *p);
struct bs *ga_mergesort_bs(struct bs *p);
struct sig *ga_mergesort_sig(struct sig *p);

#endif
