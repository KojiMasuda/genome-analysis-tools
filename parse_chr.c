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

static char *extract_val_line (const char *line, int col, const char sep);
static struct chr_block *chr_block_add (const char *chr, struct chr_block **chr_block_head);
static int chr_block_append (const char *chr, struct chr_block **chr_block_head);
static int chr_block_fa_append (struct chr_block_fa **chr_block_head, const char *chr, const char *letter);
static struct bs *bs_add (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const char strand, const char *line);
static struct sig *sig_add (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const float val);
static int ref_append (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const char strand, const char *ex_st, const char *ex_ed, const char *line);
static char *pickup_str (const char *str, const int st);

/*pointer which must be freed: char *ga_header_line */
/*
 * This is the main function for parsing.
 * *filename: input file name
 * **chr_block_head: pointer of pointer to struct chr_block.
 * col_chr: column num of chr
 * col_st: column num of start
 * col_ed: column num of end
 * col_strand: column num of strand
 * hf: header flag. If 1, header is obtained from the first line of input file and pointed by global variable, ga_header_line.
 */
void ga_parse_chr_bs (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int col_strand, int hf)
{
  char line[LINE_STR_LEN];
  char *chr, *st, *ed, *strand, *e;
  FILE *fp;
  if ((fp = fopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
//    goto err;
    exit(EXIT_FAILURE);
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
    if (line[0] == '#') continue;
    chr = extract_val_line(line, col_chr, '\t'); //extracting chr, start, end positions
    st  = extract_val_line(line, col_st , '\t');
    ed  = extract_val_line(line, col_ed , '\t');
    chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)
    if (col_strand >= 0) {
      strand = extract_val_line(line, col_strand , '\t');
      bs_add (chr, chr_block_head, strtoul(st, &e, 10), strtoul(ed, &e, 10), strand[0], line); //adding bs with strand info
      free(strand);
    } else {
      bs_add (chr, chr_block_head, strtoul(st, &e, 10), strtoul(ed, &e, 10), '.', line); //adding bs
    }
    free(chr);
    free(st);
    free(ed);
  }

  fclose(fp);
  return;

err:
  fclose(fp);
  return;
}

/*pointer which must be freed: char *ga_header_line */
/*
 * This is the main function for parsing.
 * *filename: input file name
 * **chr_block_head: pointer of pointer to struct chr_block.
 * col_chr: column num of chr
 * col_st: column num of start
 * col_ed: column num of end
 * col_strand: column num of strand
 * col_ex_st: column num of exon start
 * col_ex_ed: column num of exon end
 * hf: header flag. If 1, header is obtained from the first line of input file and pointed by global variable, ga_header_line.
 */
int ga_parse_chr_ref (const char *filename, struct chr_block **chr_block_head, int col_chr, int col_st, int col_ed, int col_strand, int col_ex_st, int col_ex_ed, int hf)
{
  char line[LINE_STR_LEN];
  char *chr, *st, *ed, *strand, *ex_st, *ex_ed, *e;
  FILE *fp;
  if ((fp = fopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
    exit(EXIT_FAILURE);
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
    if (line[0] == '#') continue;
    chr = extract_val_line(line, col_chr, '\t'); //extracting chr, start, end, ex_st, ex_ed positions
    st  = extract_val_line(line, col_st , '\t');
    ed  = extract_val_line(line, col_ed , '\t');
    ex_st  = extract_val_line(line, col_ex_st , '\t');
    ex_ed  = extract_val_line(line, col_ex_ed , '\t');
    if(chr_block_append (chr, chr_block_head) != 0){ //appending chr link list (if the chr is already linked, the input chr is just ignored)
      LOG("error: error in chr_block_append function.");
      goto err;
    }

    if (col_strand >= 0) {
      strand = extract_val_line(line, col_strand , '\t');
      if (ref_append (chr, chr_block_head, strtoul(st, &e, 10), strtoul(ed, &e, 10), strand[0], ex_st, ex_ed, line) != 0){ //appending ref
        LOG("error: error in ref_append function.");
        goto err;
      }
      free(strand);
    } else {
      if (ref_append (chr, chr_block_head, strtoul(st, &e, 10), strtoul(ed, &e, 10), '.', ex_st, ex_ed, line) != 0){ //appending ref
        LOG("error: error in ref_append function.");
        goto err;
      }
    }
    free(chr);
    free(st);
    free(ed);
    free(ex_st);
    free(ex_ed);
  }

  fclose(fp);
  return 0;

err:
  fclose(fp);
  if (chr) free(chr);
  if (st) free(st);
  if (ed) free(ed);
  if (ex_st) free(ex_st);
  if (ex_ed) free(ex_ed);
  if (strand) free(strand);
  return -1;
}

/*pointer which must be freed: char *chr, *letter*/
/*
 * This is the main function for parsing fasta file.
 * *filename: input file name
 * **chr_block_head: pointer of pointer to struct chr_block_fa.
 * *gt: genome table, pointer of struct chr_block.
 */
int ga_parse_chr_fa (const char *filename, struct chr_block_fa **chr_block_head, struct chr_block *chr_block_head_gt)
{
  char line[LINE_STR_LEN], chr_tmp[LINE_STR_LEN], *letter = NULL;
  int i, init=0;
  unsigned long j = 0; //for letter position
  struct chr_block *ch;
  FILE *fp;
  if ((fp = fopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
    exit(EXIT_FAILURE);
  }

  while (fgets(line, LINE_STR_LEN * sizeof(char), fp) != NULL) {
    if (strlen(line) >= LINE_STR_LEN -1) {
      LOG("errer: line length is too long.");
      goto err;
    }
    if (line[0] == '#') continue;
    if (line[0] == '>') { //if chromosome
      if (line[1] == '\n') {
        LOG("error: no chromosome.");
        goto err;
      }
      if (!init) init = 1; //if first chromosome
      else { //if following chromosome
        if(chr_block_fa_append(chr_block_head, chr_tmp, letter) != 0) {
          LOG("error: error in chr_block_fa_append function.");
          goto err;
        }
        free(letter);
        letter = NULL;
      }

      i = 0;
      while(line[i+1] != ' ' && line[i+1] != '\n') {
        chr_tmp[i] = line[i+1]; //copying chromosome
        i++;
      }
      chr_tmp[i] = '\0'; //null

      for (ch = chr_block_head_gt; ch; ch = ch->next) {
        if (!(strcmp(ch->chr, chr_tmp))) break;
      }
      if (ch == NULL) {
        printf("error: chromosome '%s' is not found in genome table file.\n", chr_tmp);
        goto err;
      }
      letter = (char*)calloc(ch->bs_list->st + 10000, sizeof(char)); //allocating letter for each chromosome
      j = 0; //reset j
    } else { //if letter
/*      i = 0; //this one is slow...
      while(line[i] != '\n' && line[i] != '\0') {
        let_tmp[i] = line[i]; //copying letter
        i++;
      }
      let_tmp[i] = '\0'; //null
      if (strlen(let_tmp) + strlen(letter) + 1 > ch->bs_list->st + 10000) { //if the stored letter is more than genome table
        printf("error: stored letter(%lu char) is longer than genome table data(%lu).\n", (unsigned long)(strlen(let_tmp) + strlen(letter)), ch->bs_list->st);
        goto err;
      }
      strcat(letter, let_tmp); //concatenating letters... */

/*      if (strlen(line) - 1 + strlen(letter) + 1 > ch->bs_list->st + 10000) { //if the stored letter is more than genome table (this one is also slow...)
        printf("error: stored letter(%lu char) is longer than genome table data(%lu).\n", (unsigned long)(strlen(line) - 1 + strlen(letter)), ch->bs_list->st);
        goto err;
      }
      strncat(letter, line, strlen(line) - 1); //concatenating letters... */
      i = 0;
      while(line[i] != '\n' && line[i] != '\0') {
        if (j + 2 > ch->bs_list->st + 10000) { //if the stored letter is more than genome table
          printf("error: stored letter(%lu char) is longer than genome table data(%lu).\n", j + 1, ch->bs_list->st);
          goto err;
        }
        letter[j] = line[i]; //copying letter
        i++;
        j++;
      }
    }
  } 

  if(chr_block_fa_append(chr_block_head, chr_tmp, letter) != 0) {
    LOG("error: error in chr_block_fa_append function.");
    goto err;
  }

  free(letter);
  fclose(fp);
  return 0;

err:
  if (letter) free(letter);
  fclose(fp);
  return -1;
}

/*
 * This creates random binding sites for simulation.
 * **chr_block_head   : pointer of pointer to struct chr_block.
 * *chr_block_head_ori: pointer to struct chr_block original of which summit number is used for picking up random positions.
 * *chr_table         : pointer to struct chr_block genome table to know the length of each chr.
 */
void ga_parse_chr_bs_rand (struct chr_block **chr_block_head, struct chr_block *chr_block_head_ori, struct chr_block *chr_table, int hw)
{
  unsigned long i;
  int rvalue;
  struct chr_block *ch, *c_table;
  clock_t cl;

  cl=clock();
  srand((unsigned)time(NULL)*cl); //seeds

  for (ch = chr_block_head_ori; ch; ch = ch -> next) {
    for (c_table = chr_table; c_table; c_table = c_table -> next) { //checking ch is in c_table
      if (!strcmp(ch->chr, c_table->chr)) break;
    }

    if (c_table == NULL) {
      LOG("error: chr is not in the genome table.");
      goto err;
    }

    chr_block_add (ch->chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)
    for (i=0; i < ch->bs_nb; i++) {
      rvalue = (rand()) % (c_table->bs_list->st - hw) + hw + 1; //rvalue must be 1-chr length
      if (rvalue % 2) {
        bs_add (ch->chr, chr_block_head, rvalue, rvalue, '+', "."); //adding bs
      }
      else {
        bs_add (ch->chr, chr_block_head, rvalue, rvalue, '-', ".");
      }
    }
  }

  return;

err:
  return;
}

/*
 * This is the main function for parsing.
 * *filename: input file name
 * **chr_block_head: pointer of pointer to struct chr_block.
 * col_chr: column num of chr
 * col_st: column num of start
 * col_ed: column num of end
 * hf: header flag. If 1, header is obtained from the first line of input file and pointed by global variable, ga_header_line.
 */
void ga_parse_bedgraph (const char *filename, struct chr_block **chr_block_head)
{
  char line[LINE_STR_LEN];
  char *chr, *st, *ed, *val, *e;
  FILE *fp;
  if ((fp = fopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
//    goto err;
    exit(EXIT_FAILURE);
  }

  while (fgets(line, LINE_STR_LEN * sizeof(char), fp) != NULL) {
    if (strlen(line) >= LINE_STR_LEN -1) {
      LOG("errer: line length is too long.");
      goto err;
    }
    if (line[0] == '#') continue;
    chr = extract_val_line(line, 0, '\t'); //extracting chr, start, end positions and val
    st  = extract_val_line(line, 1, '\t');
    ed  = extract_val_line(line, 2, '\t');
    val = extract_val_line(line, 3, '\t');
    chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)
    sig_add (chr, chr_block_head, strtoul(st, &e, 10), strtoul(ed, &e, 10), atof(val)); //adding bs
    free(chr);
    free(st);
    free(ed);
    free(val);
  }

  fclose(fp);
  return;

err:
  fclose(fp);
  return;
}

/*pointer which must be freed: char *strtemp */
/*
 * This extract a value from tab-delimited line
 * *line: input line
 * col  : column number for extraction
 * sep  : char that separates *line. Ex ' ' or '\t'.
 */
static char *extract_val_line (const char *line, int col, const char sep) {
  int x = 0, i = 0, j = 0;
  int len = strlen(line);
  
  char *strtemp = calloc (len, sizeof(char));
  if (strtemp == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  while (line[i] != '\0' && line[i] != '\n' && x <= col) {
    if (line[i] == sep) {
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

  /*initialization of bs and sig block*/
  p -> bs_init = 0;
  p -> sig_init = 0;
  p -> ref_init = 0;

  p -> next = *chr_block_head; //adding new chr block
  *chr_block_head = p;

  return (p);

err:
  if (p) free(p);
  return (NULL);
}

/*pointer which must be freed: struct chr_block *p, p->chr*/
/*
 * This appends new struct chr_block list
 * *chr: pointer to chromosome name
 * **chr_block_head: pointer of pointer to the head of the link
 */
static int chr_block_append (const char *chr, struct chr_block **chr_block_head)
{
  struct chr_block *p;
  struct chr_block *ch;

  for (ch = *chr_block_head; ch; ch = ch->next) {//checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) return 0;
  }

  p = malloc(sizeof(struct chr_block));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> chr = strdup(chr); //assigning chromosome name

  /*initialization of bs and sig block*/
  p -> bs_init = 0;
  p -> sig_init = 0;
  p -> ref_init = 0;

  if (*chr_block_head == NULL) { //if chromosome is the first one
    *chr_block_head = p;
  } else {
    (*chr_block_head) -> tail -> next = p; //if not parenthesis, *chr_block_head -> tail means "tail pointer" of "pointer of pointer chr_block_head", not "tail pointer" of "pointer of *chr_block_head"
  }

  (*chr_block_head) -> tail = p;
  p -> next = NULL;

  return 0;

err:
  if (p) free(p);
  return -1;
}

/*pointer which must be freed: p(if needed), p->chr, p->letter*/
/*
 * This appends new struct chr_block_fa list
 * **chr_block_head: pointer of pointer to the head of the link(structure chr_block_fa)
 * *chr: pointer to chromosome name
 * *letter: pointer to letter
 */
static int chr_block_fa_append (struct chr_block_fa **chr_block_head, const char *chr, const char *letter)
{
  struct chr_block_fa *p;

  p = malloc(sizeof(struct chr_block_fa));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> chr = strdup(chr); //assigning chromosome name
  p -> letter = strdup(letter); //assigning letter for chromosome
  p -> letter_len = (unsigned long)strlen(letter); //storing letter length
  
  if (*chr_block_head == NULL) { //if chromosome is the first one
    *chr_block_head = p;
  } else {
    (*chr_block_head) -> tail -> next = p; //if not parenthesis, *chr_block_head -> tail means "tail pointer" of "pointer of pointer chr_block_head", not "tail pointer" of "pointer of *chr_block_head"
  }

  (*chr_block_head) -> tail = p;
  p -> next = NULL;

  return 0;

err:
  if (p) free(p);
  return -1;
}

/*pointer which must be freed: struct bs *p, p->line, */
/*
 * This adds new struct bs list
 * *chr: pointer to chr name
 * **chr_block_head: pointer of pointer to the head of the link
 * st: start position
 * ed: end position
 * strand: strand either '+', '-' or '.'.
 * *line: pointer to each line which is read
 */
static struct bs *bs_add (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const char strand, const char *line)
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
  p -> strand = strand; //assigning strand info
  p -> line = strdup(line); //assigning line

  for (ch = *chr_block_head; ch; ch = ch->next) { //checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) break;
  }

  if (ch == NULL) {
    fprintf(stderr, "error: chr %s is not in the chr block list", chr);
    goto err;
  }

  /*initialization of bs*/
  if (!ch -> bs_init) { //if the bs is the first one to be added.
    ch -> bs_list = NULL;
    ch -> bs_init = 1; //initialization
    ch -> bs_nb = 1; //binding site number
  } else { //if the bs is not the first one to be added.
    ch -> bs_nb = ch -> bs_nb + 1;
  }

  p -> prev = NULL;
  p -> next = ch -> bs_list; //adding new bs
  if (ch -> bs_list) ch -> bs_list -> prev = p;
  ch -> bs_list = p;

  return (p);

err:
  if (p) free (p);
  return (NULL);
}

/*pointer which must be freed: struct bs *p, p->line, p->ex_st, p->ex_ed*/
/*
 * This appends new struct ref list
 * *chr: pointer to chr name
 * **chr_block_head: pointer of pointer to the head of the link
 * st: start position
 * ed: end position
 * strand: strand either '+', '-' or '.'.
 * *ex_st: pointer to exon start positions
 * *ex_ed: pointer to exon end positions
 * *line: pointer to each line which is read
 */
static int ref_append (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const char strand, const char *ex_st, const char *ex_ed, const char *line)
{
  struct ref *p;
  struct chr_block *ch;

  p = malloc(sizeof(struct ref));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> st = st; //assigning start position
  p -> ed = ed; //assigning end position
  p -> strand = strand; //assigning strand info
  p -> ex_st = strdup(ex_st); //assigning exon start
  p -> ex_ed = strdup(ex_ed); //assigning exon end
  p -> line = strdup(line); //assigning line
  p -> ov_gene = NULL; //at this point, ov_gene is null.

  for (ch = *chr_block_head; ch; ch = ch->next) { //checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) break;
  }

  if (ch == NULL) {
    fprintf(stderr, "error: chr %s is not in the chr block list", chr);
    goto err;
  }

  /*initialization of bs*/
  if (!(ch -> ref_init)) { //if the ref is the first one to be added.
    ch -> ref_list = p;
    p -> prev = NULL;
    ch -> ref_init = 1; //initialization
  } else { //if the bs is not the first one to be added.
    ch -> ref_list -> tail -> next = p;
    p -> prev = ch -> ref_list -> tail;
  }

  p -> next = NULL;
  ch -> ref_list -> tail = p;

  return 0;

err:
  if (p) free (p);
  return -1;
}

/*pointer which must be freed: struct sig *p */
/*
 * This adds new struct bs list
 * *chr: pointer to chr name
 * **chr_block_head: pointer of pointer to the head of the link
 * st: start position
 * ed: end position
 * *line: pointer to each line which is read
 */
static struct sig *sig_add (const char *chr, struct chr_block **chr_block_head, const unsigned long st, const unsigned long ed, const float val)
{
  struct sig *p;
  struct chr_block *ch;

  p = malloc(sizeof(struct sig));
  if (p == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  p -> st = st; //assigning start position
  p -> ed = ed; //assigning end position
  p -> val = val; //assigning value

  for (ch = *chr_block_head; ch; ch = ch->next) { //checking chr is already in chr_block list
    if (!strcmp(chr, ch->chr)) break;
  }

  if (ch == NULL) {
    fprintf(stderr, "error: chr %s is not in the chr block list", chr);
    goto err;
  }

  /*initialization of sig block*/
  if (!ch -> sig_init) { //if the sig is the first one to be added.
    ch -> sig_list = NULL;
    ch -> sig_init = 1; //initialization
  }

  p -> prev = NULL;
  p -> next = ch -> sig_list; //adding new sig
  if (ch -> sig_list) ch -> sig_list -> prev = p;
  ch -> sig_list = p;

  return (p);

err:
  if (p) free (p);
  return (NULL);
}

/*pointer which must be freed: struct sig *p */
/*
 * This parses separated wig.gz files.
 * *filename       : file name
 * **chr_block_head: pointer of pointer to struct chr_block
 */
void ga_parse_sepwiggz (const char *filename, struct chr_block **chr_block_head)
{
  char line[LINE_STR_LEN], tmpfile[128], str[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN], fileline[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN], str_last[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN];

  char *step=NULL, *chr=NULL, *chr_tmp=NULL, *span=NULL, *span_tmp=NULL, *val_tmp=NULL, *st_tmp=NULL, *start=NULL, *step_tmp=NULL, *step_p=NULL, *e;
  int span_val=1, step_val=0;
  unsigned long st=0;
  struct gzFile_s *gfp = NULL;
  FILE *fp = NULL;

  sprintf(tmpfile,"/tmp/ls%d.tmp",getpid()); //tmp file
  sprintf(str,"ls -1 %s_*.wig.gz > %s", filename, tmpfile); //list of wig.gz file is written in tmp file
  if(system(str) == -1) {
    LOG("error: system error for 'ls -1 filename*.wig.gz > tmpfile'");
    goto err;
  }

  if ((fp = fopen (tmpfile, "r")) == NULL) {
    LOG("errer: tmp file cannot be open.");
//    goto err;
    exit(EXIT_FAILURE);
  }

  while (fgets(fileline, (PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN) * sizeof(char), fp) != NULL) { //reading file name list
    if (strlen(fileline) >= PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN -1) {
      LOG("errer: filename is too long.");
      goto err;
    }

    fileline[strlen(fileline) - 1] = '\0'; // the last \n is set as \0
    if ((gfp = gzopen (fileline, "r")) == NULL) {
      LOG("errer: input file cannot be open.");
//      goto err;
      exit(EXIT_FAILURE);
    }

    while (gzgets(gfp, line, LINE_STR_LEN * sizeof(char)) != NULL) { //reading each line of each wig.gz
      if (strlen(line) >= LINE_STR_LEN -1) {
        LOG("errer: line length is too long.");
        goto err;
      }

      step = extract_val_line(line, 0, '\t'); //extracting step, chr, and so on...

      if (!strcmp(step, "variableStep")) {
        chr_tmp  = extract_val_line(line, 1, '\t');
        span_tmp = extract_val_line(line, 2, '\t');
        chr = pickup_str (chr_tmp, 6);
        chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)
        if (span_tmp!=NULL) {
          span = pickup_str (span_tmp, 5);
          span_val = atoi(span);
          if (span) free(span);
          free(span_tmp);
        }
        if (chr_tmp) free(chr_tmp);
        break;
      } else if (!strcmp(step, "fixedStep")) {
        chr_tmp  = extract_val_line(line, 1, '\t');
        st_tmp   = extract_val_line(line, 2, '\t');
        step_tmp = extract_val_line(line, 3, '\t');
        span_tmp = extract_val_line(line, 4, '\t');

        chr = pickup_str (chr_tmp, 6);
        chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)

        start = pickup_str (st_tmp, 6);
        st = strtoul(start, &e, 10); //start pos (unsigned long)
        step_p = pickup_str (step_tmp, 5);
        step_val = atoi(step_p);
        if (span_tmp!=NULL) {
          span = pickup_str (span_tmp, 5);
          span_val = atoi(span);
          if (span) free(span);
          free(span_tmp);
        }
        if (chr_tmp) free(chr_tmp);
        if (st_tmp) free(st_tmp);
        if (start) free(start);
        if (step_tmp) free(step_tmp);
        if (step_p) free(step_p);
        break;
      }
    }

    if (!strcmp(step, "variableStep")) {
      while (gzgets(gfp, line, LINE_STR_LEN * sizeof(char)) != NULL) { //reading each line
        st_tmp  = extract_val_line(line, 0, '\t');
        val_tmp = extract_val_line(line, 1, '\t');
        st = strtoul(st_tmp, &e, 10);
        sig_add (chr, chr_block_head, st, st + span_val, atof(val_tmp));
        free(st_tmp);
        free(val_tmp);
      }
    } else if (!strcmp(step, "fixedStep")) {
      while (gzgets(gfp, line, LINE_STR_LEN * sizeof(char)) != NULL) { //reading each line
        val_tmp = extract_val_line(line, 0, '\t');
        sig_add (chr, chr_block_head, st, st + span_val, atof(val_tmp));
        free(val_tmp);
        st += step_val;
      }
    }
  }

  if (chr) free(chr);
  if (step) free(step);
  gzclose(gfp);
  fclose(fp);

  sprintf(str_last,"rm -f %s",tmpfile);
  if(system(str_last) == -1) {
    LOG("error: system error for 'rm -f tmpfile'");
    goto err;
  }

  return;

err:
  gzclose(gfp);
  fclose(fp);
  if (step) free(step);
  if (chr) free(chr);
  if (chr_tmp) free(chr_tmp);
  if (span) free(span);
  if (span_tmp) free(span_tmp);
  if (val_tmp) free(val_tmp);
  if (st_tmp) free(st_tmp);
  if (start) free(start);
  if (step_tmp) free(step_tmp);
  if (step_p) free(step_p);
  return;
}

/*
 * This parses one wig.gz files.
 * *filename       : file name
 * **chr_block_head: pointer of pointer to struct chr_block
 */
void ga_parse_onewiggz (const char *filename, struct chr_block **chr_block_head)
{
  char line[LINE_STR_LEN], stephold[20];

  char *step=NULL, *chr=NULL, *chr_tmp=NULL, *span=NULL, *span_tmp=NULL, *val_tmp=NULL, *st_tmp=NULL, *start=NULL, *step_tmp=NULL, *step_p=NULL, *e;
  int span_val=1, st=0, step_val=0, fl=0;
  struct gzFile_s *gfp = NULL;


  if ((gfp = gzopen (filename, "r")) == NULL) {
    LOG("errer: input file cannot be open.");
//    goto err;
    exit(EXIT_FAILURE);
  }

  while (gzgets(gfp, line, LINE_STR_LEN * sizeof(char)) != NULL) { //reading each line
    if (strlen(line) >= LINE_STR_LEN -1) {
      LOG("errer: line length is too long.");
      goto err;
    }

    step = extract_val_line(line, 0, ' '); //extracting step, chr, and so on...

    if (!strcmp(step, "variableStep")) {
      chr_tmp  = extract_val_line(line, 1, ' ');
      span_tmp = extract_val_line(line, 2, ' ');

      chr = pickup_str (chr_tmp, 6);
      chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)

      if (span_tmp!=NULL) {
        span = pickup_str (span_tmp, 5);
        span_val = atoi(span);
        if (span) free(span);
        free(span_tmp);
      }
      if (chr_tmp) free(chr_tmp);
      fl = 1; //with the flag, the program can read the value.
      strcpy(stephold, step);
      continue;
    } else if (!strcmp(step, "fixedStep")) {
      chr_tmp  = extract_val_line(line, 1, ' ');
      st_tmp   = extract_val_line(line, 2, ' ');
      step_tmp = extract_val_line(line, 3, ' ');
      span_tmp = extract_val_line(line, 4, ' ');

      chr = pickup_str (chr_tmp, 6);
      chr_block_add (chr, chr_block_head); //adding chr link list (if the chr is already linked, the input chr is just ignored)

      start = pickup_str (st_tmp, 6);
      st = strtoul(start, &e, 10); //start pos (unsigned long)
      step_p = pickup_str (step_tmp, 5);
      step_val = atoi(step_p);
      if (span_tmp!=NULL) {
        span = pickup_str (span_tmp, 5);
        span_val = atoi(span);
        if (span) free(span);
        free(span_tmp);
      }
      if (chr_tmp) free(chr_tmp);
      if (st_tmp) free(st_tmp);
      if (start) free(start);
      if (step_tmp) free(step_tmp);
      if (step_p) free(step_p);
      fl = 1; //with the flag, the program can read the value
      strcpy(stephold, step);
      continue;
    } else if (fl && !strcmp(stephold, "variableStep")) {
      st_tmp  = extract_val_line(line, 0, ' ');
      val_tmp = extract_val_line(line, 1, ' ');
      st = strtoul(st_tmp, &e, 10);
      sig_add (chr, chr_block_head, st, st + span_val, atof(val_tmp));
      free(st_tmp);
      free(val_tmp);
    } else if (fl && !strcmp(stephold, "fixedStep")) {
      val_tmp = extract_val_line(line, 0, ' ');
      sig_add (chr, chr_block_head, st, st + span_val, atof(val_tmp));
      free(val_tmp);
      st += step_val;
    }
  }

  if (chr) free(chr);
  if (step) free(step);
  gzclose(gfp);

  return;

err:
  gzclose(gfp);
  if (step) free(step);
  if (chr) free(chr);
  if (chr_tmp) free(chr_tmp);
  if (span) free(span);
  if (span_tmp) free(span_tmp);
  if (val_tmp) free(val_tmp);
  if (st_tmp) free(st_tmp);
  if (start) free(start);
  if (step_tmp) free(step_tmp);
  if (step_p) free(step_p);
  return;
}

/*
 *pointer which must be freed: char *strtemp
 * This picks up string from str
 * *str: original string
 * st  : the number where picking starts. For example, if the str = "chr=chr1" and you want to get "chr1", st = 4. For now, we cannot pick up string like "chr1" from "chr1=chr" or "aaa=chr1,bbb".
 */
static char *pickup_str (const char *str, const int st) {
  int i;
  int len = strlen(str);

  char *strtemp = calloc (len, sizeof(char));
  if (strtemp == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  for (i = 0; ;i++) {
    strtemp[i] = str[i+st];
    if (str[i+st] == '\0') break;
  }

  return strtemp;

err:
  if (strtemp) free(strtemp);
  return (NULL);
}

/*
 * This simply sum total peak number from each chr.
 * *chr_block_head: pointer to struct chr_block
 */
unsigned long ga_count_peaks (struct chr_block *chr_block_head)
{
  unsigned long smt = 0;
  struct chr_block *ch;
  for (ch = chr_block_head; ch; ch = ch->next) {
    smt = smt + ch -> bs_nb; //adding each bs_nb for each chr
  }
  return smt;
}

/*
 * This frees struct chr_block list
 * **chr_block: pointer of pointer to the head of the link
 */
void ga_free_chr_block (struct chr_block **chr_block)
{
  struct chr_block *ch, *ch_tmp;
  struct bs *bs, *bs_tmp;
  struct sig *sig, *sig_tmp;
  struct ref *ref, *ref_tmp;

  ch = *chr_block;
  while (ch) {
    bs = ch -> bs_list;
    sig = ch -> sig_list;
    ref = ch -> ref_list;

    if (ch->bs_init) {
      while (bs) {
        free(bs->line);
        bs_tmp = bs->next;
        free(bs);
        bs = bs_tmp;
      }
    }

    if (ch->ref_init) {
      while (ref) {
        free(ref->ex_st);
        free(ref->ex_ed);
        if (NULL != ref->ov_gene) free(ref->ov_gene);
        free(ref->line);
        ref_tmp = ref->next;
        free(ref);
        ref = ref_tmp;
      }
    }

    if (ch->sig_init) {
      while (sig) {
        sig_tmp = sig->next;
        free(sig);
        sig = sig_tmp;
      }
    }

    free(ch->chr);
    ch_tmp = ch->next;
    free(ch);
    ch = ch_tmp;
  }

  return;
}

/*
 * This frees struct chr_block_fa list
 * **chr_block: pointer of pointer to the head of the link
 */
void ga_free_chr_block_fa (struct chr_block_fa **chr_block)
{
  struct chr_block_fa *ch, *ch_tmp;

  ch = *chr_block;
  while (ch) {
    free(ch->chr);
    free(ch->letter);
    ch_tmp = ch->next;
    free(ch);
    ch = ch_tmp;
  }

  return;
}


