#ifndef _GA_MY_H_
#define _GA_MY_H_

#include <stdio.h>
#include <stdlib.h>

void *my_malloc(size_t size);
void *my_realloc(void *ptr, size_t size);
void *my_calloc(size_t n, size_t size);
//void my_free(void *ptr);
int arr_delete_ul(unsigned long arr[], const int len, const int i);
int arr_insert_ul(unsigned long arr[], const int len, const int i, const unsigned long v);

#define MYFREE(p) {if(p){free(p); (p)=NULL;} }

#endif
