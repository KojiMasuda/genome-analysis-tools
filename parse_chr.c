/*
 * This program is one of the genome analysis tools.
 * This program parses peaks (binding site, bs) into chromosomes
 * Chromosome information is stored as link list.
 * BS information is stored as link list of BS which is connected with one chromosome link list
 */

#include "parse_chr.h"

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

static char *extract_val_line (const char *line, int col);
static struct chr_block *chr_block_add (const char *chr, struct chr_block **chr_block_head);
static struct bs *bs_add (const char *chr, struct chr_block **chr_block_head, const int st, const int ed, const char *line);

/*pointer which must be freed: char *ga_header_line */
/*
 * This is the main function for parsing.
 * *filename: input file name
 * **chr_block_head: pointer of pointer to struct chr_block.
 * col_chr: column num of chr
 * col_st: column num of start
 * col_ed: column num of end
 * hf: header flag. If 1, header is obtained from the first line of input file and pointed by global variable, ga_header_line.
 */
void ga_parse_chr_bs (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int hf)
{
  char line[LINE_STR_LEN];
  char *chr, *st, *ed;
  FILE *fp;
  if ((fp = fopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
    goto err;
  }

  if (hf) {
    if (fgets(line, LINE_STR_LEN * sizeof(char), fp) != NULL) {
      if (strlen(line) >= LINE_STR_LEN -1) {
        LOG("errer: line length is too long.");
        goto err;
      }
      else if (ga_header_line == NULL){ //so, header is obtained from a file which is opened FIRST during the whole program, basically
        ga_header_line = strdup(line);
      }
    }
  }

  while (fgets(line, LINE_STR_LEN * sizeof(char), fp) != NULL) {
    if (strlen(line) >= LINE_STR_LEN -1) {
      LOG("errer: line length is too long.");
      goto err;
    }
    chr = extract_val_line(line, col_chr); //extracting chr, start, end positions
    st = extract_val_line(line, col_st);
    ed = extract_val_line(line, col_ed);
    chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)
    bs_add (chr, chr_block_head, atoi(st), atoi(ed), line); //adding bs
    free(chr);
    free(st);
    free(ed);
  }

  fclose(fp);
  return;

err:
  return;
}

/*pointer which must be freed: char *strtemp */
/*
 * This extract a value from tab-delimited line
 * *line: input line
 * col: column number for extraction
 */
static char *extract_val_line (const char *line, int col) {
  int x = 0, i = 0, j = 0;
  int len = strlen(line);
  
  char *strtemp = calloc (len, sizeof(char));
  if (strtemp == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  while (line[i] != '\0' && line[i] != '\n' && x <= col) {
    if (line[i] == '\t') {
      i++;
      x++;
    }
    else if (col == x) { //copying the value
      strtemp[j++] = line[i++];
    }
    else {
      i++;
    }
  }
  strtemp[j] = '\0';

  if (!strlen(strtemp)) {
    LOG("error: column number is out of range.");
    goto err;
  }

  return (strtemp);

err:
  if (strtemp) free(strtemp);
  return (NULL);
}

/*pointer which must be freed: struct chr_block *p, p->chr*/
/*
 * This adds new struct chr_block list
 * *chr: pointer to chromosome name
 * **chr_block_head: pointer of pointer to the head of the link
 */
static struct chr_block *chr_block_add (const char *chr, struct chr_block **chr_block_head)
{
  struct chr_block *p;
  struct chr_block *ch;

  for (ch = *chr_block_head; ch; ch = ch->next) {//checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) return (NULL);
  }

  p = malloc(sizeof(struct chr_block));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> chr = strdup(chr); //assigning chromosome name

  p -> next = *chr_block_head; //adding new chr block
  *chr_block_head = p;

  return (p);

err:
  if (p) free(p);
  return (NULL);
}

/*pointer which must be freed: struct bs *p, p->line, */
/*
 * This adds new struct bs list
 * *chr: pointer to chr name
 * **chr_block_head: pointer of pointer to the head of the link
 * st: start position
 * ed: end position
 * *line: pointer to each line which is read
 */
static struct bs *bs_add (const char *chr, struct chr_block **chr_block_head, const int st, const int ed, const char *line)
{
  struct bs *p;
  struct chr_block *ch;

  p = malloc(sizeof(struct bs));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> st = st; //assigning start position
  p -> ed = ed; //assigning end position
  p -> line = strdup(line); //assigning line

  for (ch = *chr_block_head; ch; ch = ch->next) { //checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) break;
  }

  if (ch == NULL) {
    fprintf(stderr, "error: chr %s is not in the chr block list", chr);
    goto err;
  }

  p -> next = ch -> bs_list; //adding new bs
  ch -> bs_list = p;

  return (p);

err:
  if (p) free (p);
  return (NULL);
}


/*
 * This frees struct chr_block list
 * **chr_block: pointer of pointer to the head of the link
 */
void ga_free_chr_block (struct chr_block **chr_block)
{
  struct chr_block *ch, *ch_tmp;
  struct bs *bs, *bs_tmp;

  ch = *chr_block;
  while (ch) {
    bs = ch -> bs_list;
    while (bs) {
      free(bs->line);
      bs_tmp = bs->next;
      free(bs);
      bs = bs_tmp;
    }
    free(ch->chr);
    ch_tmp = ch->next;
    free(ch);
    ch = ch_tmp;
  }

  return;
}


