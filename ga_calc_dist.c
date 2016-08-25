/*
 * This program is one of the genome analysis tools.
 * This program calculates the distance between two peak sets or inter-summit distance (not supported yet...).
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "argument.h"
#include "sort_list.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

static void cmp_two_peaks (struct bs *bs1, struct bs *bs2, struct output **output_head);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_calc_dist\n\
Summary: report the distance of two peaks(mode:two) or inter-summit distance(mode:isd)\n\n\
Usage:   ga_calc_dist [options] -1 <file1> -2 <file2> -mode two\n\
  or:    ga_calc_dist [options] -1 <file1> -mode isd\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the header of file1 is preserved.\n\
         --col_chr1 <int>: column number for chromosome of file1 (default:0).\n\
         --col_chr2 <int>: column number for chromosome of file2 (default:0).\n\
         --col_smt1 <int>: column number for peak summit position of file1 (default:3).\n\
         --col_smt2 <int>: column number for peak summit position of file2 (default:3).\n");
  exit(0);
}

static void version()
{
  printf("'ga_calc_dist' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}

static int hf = 0;
static char *file1 = NULL;
static char *file2 = NULL;
static char *mode = NULL;
static int col_chr1 = 0, col_chr2 = 0;
static int col_smt1 = 3, col_smt2 = 3;
char *ga_header_line = NULL; //header line. Note this is global variable
static char ga_line_out[LINE_STR_LEN]; //output line with overlapping flag


static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage  },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage  },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version},
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf  },
  {"-1"          , ARGUMENT_TYPE_STRING  , &file1},
  {"-2"          , ARGUMENT_TYPE_STRING  , &file2},
  {"-mode"          , ARGUMENT_TYPE_STRING  , &mode},
  {"--col_chr1"  , ARGUMENT_TYPE_INTEGER , &col_chr1   },
  {"--col_chr2"  , ARGUMENT_TYPE_INTEGER , &col_chr2   },
  {"--col_smt1", ARGUMENT_TYPE_INTEGER , &col_smt1   },
  {"--col_smt2", ARGUMENT_TYPE_INTEGER , &col_smt2   },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL   },
};

int main (int argc, char *argv[])
{
  struct chr_block *chr_block_head1 = NULL; //for peak1
  struct chr_block *chr_block_head2 = NULL; //for peak2

  struct output *output_head = NULL; //for output

  struct chr_block *ch1, *ch2; //for "for loop of chr"
  struct bs *bs_nonov, *bs1; //for bs of file1 of which chr is not in file2 and for loop

  /*path, filename, and extension*/
  char path1[PATH_STR_LEN];
  char fn1[FILE_STR_LEN];
  char ext1[EXT_STR_LEN];
  char path2[PATH_STR_LEN];
  char fn2[FILE_STR_LEN];
  char ext2[EXT_STR_LEN];
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN]; //output file name
  time_t timer;

  argument_read(&argc, argv, args);//reading arguments
  if (file1 == NULL || mode == NULL) usage();

  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
Input file1:                     %s\n\
Input file2:                     %s\n\
mode:                            %s\n\
File1 column of chr, summit: %d, %d\n\
File2 column of chr, summit: %d, %d\n\
header flag:                     %d\n\
time:                            %s\n",\
 "ga_calc_dist", file1, file2, mode, col_chr1, col_smt1, col_chr2, col_smt2, hf, ctime(&timer) );

  if (!strcmp(mode, "two")) { //mode == two
    ga_parse_chr_bs(file1, &chr_block_head1, col_chr1, col_smt1, col_smt1, -1, hf); //parsing each binding sites for each chromosome without strand info
    ga_parse_chr_bs(file2, &chr_block_head2, col_chr2, col_smt2, col_smt2, -1, hf);

    chr_block_head1 = ga_mergesort_chr(chr_block_head1); //sorting chr block
    chr_block_head2 = ga_mergesort_chr(chr_block_head2); //sorting chr block

    for (ch1 = chr_block_head1; ch1; ch1 = ch1 -> next) { //sorting bs
      ch1 -> bs_list = ga_mergesort_bs(ch1 -> bs_list);
    }
    for (ch2 = chr_block_head2; ch2; ch2 = ch2 -> next) { //sorting bs
      ch2 -> bs_list = ga_mergesort_bs(ch2 -> bs_list);
    }


    for (ch1 = chr_block_head1; ch1; ch1 = ch1->next) {
      for (ch2 = chr_block_head2; ch2; ch2 = ch2->next) {
        if (!strcmp(ch1->chr, ch2->chr)) { //comparing peaks on the same chr
          cmp_two_peaks (ch1->bs_list, ch2->bs_list, &output_head);
          break;
        }
      }
      if (ch2 == NULL) { //if chr_block2 doesn't have ch1
        for (bs_nonov = ch1->bs_list; bs_nonov; bs_nonov = bs_nonov->next) {
          if (add_one_val(ga_line_out, bs_nonov->line, "NA\tNA\n") < 0){
            LOG("error: output line was too long.");
            goto err; //making output link list with NA
          }
          ga_output_append(&output_head, ga_line_out);
        }
      }
    }

    ga_parse_file_path (file1, path1, fn1, ext1); //parsing input file name into path, file name, and extension
    ga_parse_file_path (file2, path2, fn2, ext2);

    sprintf(output_name, "%s%s_closest_%s%s", path1, fn1, fn2, ext1);
    if (ga_header_line != NULL) { //if header line
      if (add_one_val(ga_line_out, ga_header_line, "dist\tsummit2\n") < 0){
        LOG("error: output line was too long.");
        goto err; //adding one extra column
      }
      ga_write_lines (output_name, output_head, ga_line_out); //note that header is line_out, not ga_header_line
    }
    else ga_write_lines (output_name, output_head, ga_header_line);

    if (ga_header_line) free(ga_header_line);
    ga_free_chr_block(&chr_block_head1);
    ga_free_chr_block(&chr_block_head2);

    ga_free_output(&output_head);

    return 0;
  } else if (!strcmp(mode, "isd")) { //mode == two or mode == isd
    ga_parse_chr_bs(file1, &chr_block_head1, col_chr1, col_smt1, col_smt1, -1, hf); //parsing each binding sites for each chromosome without strand info

    chr_block_head1 = ga_mergesort_chr(chr_block_head1); //sorting chr block

    for (ch1 = chr_block_head1; ch1; ch1 = ch1 -> next) { //sorting bs
      ch1 -> bs_list = ga_mergesort_bs(ch1 -> bs_list);
    }

    for (ch1 = chr_block_head1; ch1; ch1 = ch1->next) {
      for (bs1 = ch1->bs_list; bs1; bs1 = bs1->next) {
        if (bs1->next == NULL) break; //if the bs is the tail one
        sprintf(ga_line_out, "%s\t%lu\n", ch1->chr, bs1->next->st - bs1->st);
        ga_output_append(&output_head, ga_line_out);
      }
    }

    ga_parse_file_path (file1, path1, fn1, ext1); //parsing input file name into path, file name, and extension

    sprintf(output_name, "%s%s_ISD%s", path1, fn1, ext1);
    ga_write_lines (output_name, output_head, "chromosome\tInter_Summit_Distance\n");

    if (ga_header_line) free(ga_header_line);
    ga_free_chr_block(&chr_block_head1);
    ga_free_output(&output_head);

    return 0;
  } else { //mode == two or isd
    printf("error: input correct mode 'two' or 'isd'. Your mode was: %s.\n", mode);
    goto err;
  }

err:
  if (ga_header_line) free(ga_header_line);
  if (chr_block_head1) ga_free_chr_block(&chr_block_head1);
  if (chr_block_head2) ga_free_chr_block(&chr_block_head2);

  if (output_head) ga_free_output(&output_head);
  return -1;
}

/*
 * This function compares overlapping of bs.
 * *bs1: pointer to bs link list for file1
 * *bs2: pointer to bs link list for file2
 * **output_head: pointer of pointer to output link list for output file with flag
 */
static void cmp_two_peaks (struct bs *bs1, struct bs *bs2, struct output **output_head)
{
  struct bs *i, *j, *j_tmp;
  char tmp[64];

  for (i = bs1; i; i = i->next) {
    if (i->st <= bs2->st) { //if bs1 is the left side among all bs2
      sprintf(tmp, "%lu\t%lu\n", bs2->st - i->st, bs2->st); //distance and peak summit2
      if (add_one_val(ga_line_out, i->line, tmp) < 0) {
        LOG("error: output line was too long.");
        goto err;
      }
      ga_output_append(output_head, ga_line_out);
      continue;
    }
    for (j = bs2; j; j = j->next) {
      if (i->st < j->st) { //if bs2 exceeds bs1
        if (j->st - i->st < i->st - j->prev->st) sprintf(tmp, "%lu\t%lu\n", j->st - i->st, j->st); //if bs1 is close to j->st than (j-1)->st
        else sprintf(tmp, "%lu\t%lu\n", i->st - j->prev->st, j->prev->st);
        if (add_one_val(ga_line_out, i->line, tmp) < 0){
          LOG("error: output line was too long.");
          goto err;
        }
        ga_output_append(output_head, ga_line_out);
        break;
      }
      j_tmp = j;
    }
    if (j == NULL) { //if i is the right side among all bs2
      sprintf(tmp, "%lu\t%lu\n", i->st - j_tmp->st, j_tmp->st); //distance and peak summit2
      if (add_one_val(ga_line_out, i->line, tmp) < 0){
        goto err;
        LOG("error: output line was too long.");
      }
      ga_output_append(output_head, ga_line_out);
    }
  }
  return;

err:
  return;
}

