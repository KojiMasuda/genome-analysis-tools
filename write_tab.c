/*
 * This program is one of the genome analysis tools.
 * The job of the program is to write output lines.
 * This also get path, filename, extension from input file name
 * so that set output file name.
 */

#include "write_tab.h"
#include "parse_chr.h"
#include "ga_my.h"

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

static void get_path (char *str, const char *delim, char *path, size_t path_len, char *fn, size_t fn_len);

/*pointer which must be freed: struct output *p, p->line */
/*
 * This adds the new struct output list(bottom is head, top is tail direction...)
 * **out_head: pointer of pointer to struct ouput head list
 * *line: pointer to output line
 */
struct output *ga_output_add (struct output **out_head, const char *line)
{
  struct output *p;
  
  p = my_malloc(sizeof(struct output));
  p -> line = strdup(line);

  p -> next = *out_head;
  *out_head = p;

  return (p);
}

/*pointer which must be freed: struct output *p, p->line */
/*
 * This appends the new struct output list(bottom is tail, top is head direction!)
 * **out_head: pointer of pointer to struct ouput head list
 * *line: pointer to output line
 */
void ga_output_append (struct output **out_head, const char *line)
{
  struct output *p;
  
  p = my_malloc(sizeof(struct output));
  p -> line = strdup(line);

  if (*out_head == NULL) { //if the line is the first one
    *out_head = p;
  } else {
//    (**out_head).tail -> next = p; //can be
    (*out_head) -> tail -> next = p; //if not parenthesis, *out_head -> tail means "tail pointer" of "pointer of pointer out_head", not "tail pointer" of "pointer of *outhead"
  }

//  (**out_head).tail = p; //can be
  (*out_head) -> tail = p;
  p -> next = NULL;

  return;
}

/*
 * This parses input filename with path into path("/aaa/bbb/ccc/"), filename("ddd.eee"), extension(".txt").
 * *file_path: pointer to filename with path
 * *path: pointer to path.
 * *fn: pointer to filename.
 * *ext: pointer to extension.
 */
void ga_parse_file_path (char *file_path, char *path, char *fn, char *ext)
{
  char tmp[FILE_STR_LEN] = {0};
  char prev = '.', next;
  int i;

  memset(path, '\0', PATH_STR_LEN * sizeof(char));
  memset(fn, '\0', FILE_STR_LEN * sizeof(char));
  memset(ext, '\0', EXT_STR_LEN * sizeof(char));
  memset(tmp, '\0', FILE_STR_LEN * sizeof(char));

  get_path (file_path, "/", path, PATH_STR_LEN, tmp, FILE_STR_LEN); //getting path
  if (strchr(tmp, '.') == NULL) { //if . is not found in tmp(file name)
    strcpy(fn, tmp);
    return; //no ext
  } else {
    get_path (tmp, ".", fn, FILE_STR_LEN, ext, EXT_STR_LEN); //getting filename and extension
  }

  fn[strlen(fn) - 1] = '\0'; //deleting the last '.'
  for (i = 0; i < strlen(ext); i++) { //converting "txt" into ".txt"
    next = ext[i];
    ext[i] = prev;
    prev = next;
  }
  ext[i] = prev;

  return;
}

/*pointer which must be freed: strp, tmp*/
/*
 * Path and filename are obtained by this function. If *str is "/aaa/bbb/ccc/ddd.eee.txt", *delim="/", *path is "/aaa/bbb/ccc/", and *fn is "ddd.eee.txt"
 * If *str is "ddd.eee.txt", *delim=".", *path is "ddd.eee.", and *fn is "txt"
 * *str: pointer to input filename
 * *delim: pointer to deliminator
 * *path: pointer to path
 * path_len: MAX length for *path
 * *fn: pointer to filename
 * fn_len: MAX length for *fn
 */
static void get_path (char *str, const char *delim, char *path, size_t path_len, char *fn, size_t fn_len)
{
  char *strp = NULL;
  char *token, *tmp = NULL;

  strp = strdup(str);

  token = strsep(&strp, delim);
  while (token != NULL) {
    tmp = strdup(token); //copying previous token
    token = strsep(&strp, delim);

    if (token == NULL) break; //if str cannot be separated by delim anymore
    else if (strlen(path) + strlen(tmp) + 2 < path_len * sizeof(char)) { //concatenating path
      strncat(path, tmp, strlen(tmp) * sizeof(char));
      strncat(path, delim, sizeof(char));
    }
    else {
      fprintf(stderr, "path string length is over %lu.\n", path_len * sizeof(char));
      goto err;
    }
    my_free(tmp);
    tmp = NULL;
  }

  if (strlen(tmp) + 1 < fn_len * sizeof(char)) { //if token reaches NULL, tmp is the last token which is the filename
    strncpy(fn, tmp, strlen(tmp) * sizeof(char));
  }
  else {
    fprintf(stderr, "filename string length is over %lu.\n", fn_len * sizeof(char));
    goto err;
  }

  my_free(tmp);
  my_free(strp);

  return;

err:
  if (strp) my_free(strp);
  if (tmp) my_free(tmp);

  return;
}

/*
 * This writes the output lines.
 * *output: pointer to output filename
 * *out_head: pointer to struct ouput head link
 * *header: pointer to header
 */
void ga_write_lines (const char *output, struct output *out_head, const char *header)
{
  FILE *fp;
  struct output* o;

  if ((fp = fopen(output, "w")) == NULL) {
    LOG("error: output file cannot be open.");
    goto err;
  }

  if (header != NULL) {
    if (fputs(header, fp) == EOF) {
      LOG("error: file writing error.");
      goto err;
    }
  }

  for (o = out_head; o; o = o->next) {
    if (fputs(o->line, fp) == EOF) {
      LOG("error: file writing error.");
      goto err;
    }
  }

  fclose(fp);

err:
  return;
}

/*
 * This function add one more value to string. If line_out[xxx], line = "aaa\tbbb\n", val = "ccc\n", line_out is "aaa\tbbb\tccc\n".
 * line_out[]: char array. This must have size of LINE_STR_LEN.
 * *line: pointer to char to be added.
 * *val: pointer to char for adding. Put '\n' at the last position if you need.
 */
int add_one_val (char line_out[], const char *line, const char *val)
{
  sprintf(line_out, "%s", line);
  line_out[strlen(line_out) - 1] = '\t';
  if (strlen(line_out) + strlen(val) + 1 < LINE_STR_LEN) strncat(line_out, val, strlen(val));
  else {
    LOG("error: the output line length is too long.");
    return -1;
  }

  return 0;
}

/*
 * This frees struct ouput link list
 * **out_head: pointer of pointer to struct output list
 */
void ga_free_output (struct output **out_head)
{
  struct output *o, *tmp;

  o = *out_head;
  while (o) {
    my_free(o->line);
    tmp = o->next;
    my_free(o);
    o = tmp;
  }

  return;
}
