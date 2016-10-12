/*
 * This program is one of the genome analysis tools.
 * This program calculates the nucleotide composition in the specified reigons.
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "argument.h"
#include "sort_list.h"
#include "ga_my.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)


static void usage()
{
  printf("Tool:    ga_nuc_region\n\n\
Summary: report the nucleotide composition(number of ATCG) and AT content.\n\n\
Usage:   ga_nuc_region -fa <fasta> -region <region file> -gt <genome table>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the header of region file is preserved (default:off).\n\
         --col_chr: <int> column number for chromosome of region file (default:0).\n\
         --col_start: <int> column number for start position of region file (default:1).\n\
         --col_end: <int> column number for end position of region file (default:2).\n");
  exit(0);
}

static void version()
{
  printf("'ga_nuc_region' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}

char *ga_header_line = NULL; //header line. Note this is external global variable
static char *rgn = NULL; //region file
static char *fa = NULL; //fasta
static char *gt = NULL; //genome table
static char ga_line_out[LINE_STR_LEN] = {0}; //output line with nucleotide composition
static int col_chr = 0;
static int col_st = 1;
static int col_ed = 2;
static int hf = 0;
static char hfs[4] = "off\0";

static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage   },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage   },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version },
  {"-region"         , ARGUMENT_TYPE_STRING  , &rgn},
  {"-fa"         , ARGUMENT_TYPE_STRING  , &fa     },
  {"-gt"         , ARGUMENT_TYPE_STRING  , &gt     },
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf     },
  {"--col_chr"   , ARGUMENT_TYPE_INTEGER , &col_chr},
  {"--col_start" , ARGUMENT_TYPE_INTEGER , &col_st },
  {"--col_end"   , ARGUMENT_TYPE_INTEGER , &col_ed },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL    },
};

int main (int argc, char *argv[])
{
  /*path, filename, and extension*/
  char path[PATH_STR_LEN] = {0};
  char fn[FILE_STR_LEN] = {0};
  char ext[EXT_STR_LEN] = {0};
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN] = {0}; //output file name
  char tmp[LINE_STR_LEN] = {0};
  char *frag; //fragment

  struct chr_block *chr_block_head_rgn = NULL; //for parsing region file
  struct chr_block *ch1 ; //for 
  struct chr_block *chr_block_head_gt = NULL; //for parsing genome table
  struct chr_block_fa *chr_block_head_fa = NULL; //for parsing fasta
  struct chr_block_fa *ch2; //for
  struct bs *bs; //for
  struct output *out_head; //output header

  unsigned long i, c_A, c_T, c_G, c_C, max_len = 1; //count of ATCG and max len

  time_t timer;

  argument_read(&argc, argv, args);//reading arguments
  if (fa == NULL || rgn == NULL || gt == NULL) usage();

  if(hf) strcpy(hfs, "on\0");
  time(&timer);
  printf("Tool:                            %s\n\n\
fasta file:                      %s\n\
region file:                     %s\n\
genome table:                    %s\n\
header?:                         %s\n\
col of chr, start, end:          %d, %d, %d\n\
time:                            %s\n",\
 "ga_nuc_region", fa, rgn, gt, hfs, col_chr, col_st, col_ed, ctime(&timer) );

  ga_parse_chr_bs (rgn, &chr_block_head_rgn, col_chr, col_st, col_ed, -1, hf); //parsing region file
  // sorting summit and sig
  chr_block_head_rgn = ga_mergesort_chr(chr_block_head_rgn);

  for (ch1 = chr_block_head_rgn; ch1; ch1 = ch1 -> next) {
    ch1 -> bs_list = ga_mergesort_bs(ch1 -> bs_list);
    for (bs = ch1 -> bs_list; bs; bs = bs -> next) { //finding the max region length
      if (max_len < bs->ed - bs->st) max_len = bs->ed - bs->st;
    }
  }
  printf("\nmaximum region length: %lu\n", max_len);

  ga_parse_chr_bs (gt, &chr_block_head_gt, 0, 1, 1, -1, 0); //parsing genome tible
  if (ga_parse_chr_fa(fa, &chr_block_head_fa, chr_block_head_gt) != 0){
    LOG("error: error in ga_parse_chr_fa function.");
    goto err;
  }
  ga_parse_file_path (rgn, path, fn, ext); //parsing input file name into path, file name, and extension

  frag = (char *)my_malloc(sizeof(char) * (max_len + 2)); //fragment DNA from genome
  for (ch1 = chr_block_head_rgn; ch1; ch1 = ch1->next) {
    for (ch2 = chr_block_head_fa; ch2; ch2 = ch2->next) {
      if (!(strcmp(ch1->chr, ch2->chr))) break;
    }

    if (ch2 == NULL) {
      printf("warning: chromosome %s was not found in the fasta file, skipped\n", ch1->chr);
      continue;
    }
    else printf("calculating on %s \n", ch1->chr);

    for (bs = ch1 -> bs_list; bs; bs = bs -> next) { //fetching the nucleotid of the region
      memset(frag, '\0', (max_len + 2) * sizeof(char)); //this may slow down the program, but safe I believe...
      strncpy(frag, ch2->letter + bs->st - 1, bs->ed - bs->st + 1);  //sliced fragment DNA from genome. Starting from letter[0], win nucleotide is copied to frag. Next the slice starts from letter[0 + bs->st - 1].
      frag[bs->ed - bs->st + 1] = '\0'; //null
      if (strchr(frag, 'N') != NULL) printf("warning: N nucleotide was found in region %lu-%lu on %s. N was not counted.\n", bs->st, bs->ed, ch1->chr); //if N is found
      i = 0;
      c_A = 0;
      c_T = 0;
      c_C = 0;
      c_G = 0;
      while (frag[i] != '\0') {
        if (frag[i] == 'A') c_A++;
        else if (frag[i] == 'T') c_T++;
        else if (frag[i] == 'C') c_C++;
        else if (frag[i] == 'G') c_G++;
        else if (frag[i] == 'N') {
          i++;
          continue;
        }
        else printf("Oops! Unexpected nucleotide %c was found in %lu-%lu on %s.\n", frag[i], bs->st, bs->ed, ch1->chr);
        i++;
      }
      sprintf(tmp, "%lu\t%lu\t%lu\t%lu\t%.3lf\n", c_A, c_T, c_C, c_G, (double)(c_A + c_T)/(double)(c_A + c_T + c_C + c_G));
      if (add_one_val (ga_line_out, bs->line, tmp) != 0) {
        LOG("error: error in add_one_val function.");
        goto err;
      }
      ga_output_append (&out_head, ga_line_out);
    } //bs
  } //chromosome

  sprintf(output_name, "%s%s_nuc.txt", path, fn);

  if (ga_header_line != NULL) { //if header
    if (add_one_val (ga_line_out, ga_header_line, "count_A\tcount_T\tcount_C\tcount_G\tAT_content\n") != 0) {
      LOG("error: error in add_one_val function.");
      goto err;
    }
    ga_write_lines (output_name, out_head, ga_line_out);
  }
  else ga_write_lines (output_name, out_head, ga_header_line);

  MYFREE (ga_header_line);
  ga_free_chr_block_fa(&chr_block_head_fa);
  ga_free_chr_block(&chr_block_head_rgn);
  ga_free_chr_block(&chr_block_head_gt);
  MYFREE (frag);
  return 0;

err:
  MYFREE (ga_header_line);
  if (chr_block_head_fa) ga_free_chr_block_fa(&chr_block_head_fa);
  if (chr_block_head_rgn) ga_free_chr_block(&chr_block_head_rgn);
  if (chr_block_head_gt) ga_free_chr_block(&chr_block_head_gt);
  MYFREE (frag);
  return -1;
}

