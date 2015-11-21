#include "sort_list.h"

static struct chr_block *merge_chr(struct chr_block *a, struct chr_block *b);
static struct bs *merge_bs(struct bs *a, struct bs *b);
static struct sig *merge_sig(struct sig *a, struct sig *b);

struct chr_block *ga_mergesort_chr(struct chr_block *p)
{
  struct chr_block *a, *b;

  if (p == NULL) return p;
  if (p->next == NULL) return p;

  a = p;
  b = p->next->next;

  while (b != NULL) {
    a = a->next;
    b = b->next;
    if (b != NULL) b = b->next;
  }

  b = a->next;
  a->next = NULL;

  return merge_chr(ga_mergesort_chr(p),ga_mergesort_chr(b));
}

struct bs *ga_mergesort_bs(struct bs *p)
{
  struct bs *a, *b;

  if (p == NULL) return p;
  if (p->next == NULL) return p;

  a = p;
  b = p->next->next;

  while (b != NULL) {
    a = a->next;
    b = b->next;
    if (b != NULL) b = b->next;
  }

  b = a->next;
  a->next = NULL;

  return merge_bs(ga_mergesort_bs(p),ga_mergesort_bs(b));
}

struct sig *ga_mergesort_sig(struct sig *p)
{
  struct sig *a, *b;

  if (p == NULL) return p;
  if (p->next == NULL) return p;

  a = p;
  b = p->next->next;

  while (b != NULL) {
    a = a->next;
    b = b->next;
    if (b != NULL) b = b->next;
  }

  b = a->next;
  a->next = NULL;

  return merge_sig(ga_mergesort_sig(p),ga_mergesort_sig(b));
}

static struct chr_block *merge_chr(struct chr_block *a, struct chr_block *b)
{
  struct chr_block *x, head;

  x = &head;

  while (a && b) {
      if (strcmp(a->chr,b->chr)<0) {
          x->next = a;
          a = a->next;
          x = x->next;
      }
      else {
          x->next = b;
          b = b->next;
          x = x->next;
      }
  }

  if (a == NULL)
    x->next = b;
  else
    x->next = a;

  return head.next;
}


static struct bs *merge_bs(struct bs *a, struct bs *b)
{
  struct bs *x, head;

  x = &head;

  while (a && b) {
      if (a->st <= b->st) {
          x->next = a;
          a = a->next;
          x = x->next;
      }
      else {
          x->next = b;
          b = b->next;
          x = x->next;
      }
  }

  if (a == NULL)
    x->next = b;
  else
    x->next = a;

  return head.next;
}

static struct sig *merge_sig(struct sig *a, struct sig *b)
{
  struct sig *x, head;

  x = &head;

  while (a && b) {
      if (a->st <= b->st) {
          x->next = a;
          a = a->next;
          x = x->next;
      }
      else {
          x->next = b;
          b = b->next;
          x = x->next;
      }
  }

  if (a == NULL)
    x->next = b;
  else
    x->next = a;

  return head.next;
}

