/*
 * This is a series of my functions.
 */
#include "ga_my.h"

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)


void *my_malloc(size_t size)
{
  void *p = NULL;
  p = malloc(size);
  if (NULL == p) {
    LOG("error: failed malloc.");
    exit(EXIT_FAILURE);
  }

  return p;
}

void *my_realloc(void *ptr, size_t size)
{
  void *p = NULL;
  p = realloc(ptr, size);
  if (NULL == p) {
    LOG("error: failed realloc.");
    free(ptr);
    exit(EXIT_FAILURE);
  }

  return p;
}

void *my_calloc(size_t n, size_t size)
{
  void *p;

  p = calloc(n, size);
  if (NULL == p) {
    LOG("error: failed calloc.");
    exit(EXIT_FAILURE);
  }

  return p;
}

/*
 * this deletes element i of unsigned long arr
 * arr: array where element i is deleted
 * len: length of the array. Note that user must check the len is less than or equal to arr length
 * i: element number
 */
int arr_delete_ul(unsigned long arr[], const int len, const int i)
{
  int x;

  for (x = i; x < len - 1; x++) {
    arr[x] = arr[x + 1];
  }
  arr[x] = 0; //reset the end element of the arr

  return 0;
}

/*
 * this inserts element i of unsigned long arr
 * arr: array where v is inserted to element i
 * len: length of the original array. Note that user must check the len + 1 is less than or equal to arr length
 * i: element number
 * v: element value
 */
int arr_insert_ul(unsigned long arr[], const int len, const int i, const unsigned long v)
{
  int x;

  for (x = len; x > i; x--) {
    arr[x] = arr[x-1];
  }
  arr[x] = v; //inserting v

  return 0;
}














