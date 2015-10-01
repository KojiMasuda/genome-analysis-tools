#ifndef _WRITE_TAB_H_
#define _WRITE_TAB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PATH_STR_LEN 1024 //path length
#define FILE_STR_LEN 256 //filename length
#define EXT_STR_LEN 80 //extension length

/*
 * struct output
 * This is a link list
 * each line for output is linked
 */
struct output {
  char *line;
  struct output *next;
};

struct output *ga_output_add (struct output **out_head, const char *line);
void ga_free_output (struct output **out_head);
void ga_parse_file_path (char *file_path, char *pathp, char *fnp, char *extp);
void ga_write_lines (const char *output, struct output *out_head, const char *header);

#endif
