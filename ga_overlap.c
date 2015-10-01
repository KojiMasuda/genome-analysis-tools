/*
 * This program is one of the genome analysis tools.
 * This program check the overlapping of two peaks/summits.
 * The three files are saved. File1_over_File2.txt, File1_nonover_File2.txt, and File1_vs_File2.txt.
 * In File1_vs_File2.txt, flags for OV/NONOV are added at the last+1 column.
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "argument.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

static void cmp_overlap (struct bs *bs1, struct bs *bs2, struct output **output_head, struct output **ov_head, struct output **nonov_head);
static int add_one_val (char line_out[], const char *line, const char *val);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Summary: report the overlaps of two peaks/summits\n\n\
Usage:   ga_overlap [options] -1 <file1> -2 <file2>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the header of file1 is preserved.\n\
         --col_chr1 <int>: column number for chromosome of file1 (default:0).\n\
         --col_chr2 <int>: column number for chromosome of file2 (default:0).\n\
         --col_start1 <int>: column number for peak start position of file1 (default:1).\n\
         --col_start2 <int>: column number for peak start position of file2 (default:1).\n\
         --col_end1 <int>: column number for peak end position of file1 (default:2).\n\
         --col_end2 <int>: column number for peak end position of file2 (default:2).\n");
  exit(0);
}

static void version()
{
  printf("'ga_overlap' in genome analysis tools version: 0.0.1\n");
  exit(0);
}

static int hf = 0;
static char *file1 = NULL;
static char *file2 = NULL;
static int col_chr1 = 0, col_chr2 = 0;
static int col_st1 = 1, col_st2 = 1;
static int col_ed1 = 2, col_ed2 = 2;
char *ga_header_line = NULL; //header line. Note this is global variable
static char ga_line_out[LINE_STR_LEN]; //output line with overlapping flag

static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage  },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage  },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version},
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf  },
  {"-1"          , ARGUMENT_TYPE_STRING  , &file1},
  {"-2"          , ARGUMENT_TYPE_STRING  , &file2},
  {"--col_chr1"  , ARGUMENT_TYPE_INTEGER , &col_chr1   },
  {"--col_chr2"  , ARGUMENT_TYPE_INTEGER , &col_chr2   },
  {"--col_start1", ARGUMENT_TYPE_INTEGER , &col_st1   },
  {"--col_start2", ARGUMENT_TYPE_INTEGER , &col_st2   },
  {"--col_end1"  , ARGUMENT_TYPE_INTEGER , &col_ed1   },
  {"--col_end2"  , ARGUMENT_TYPE_INTEGER , &col_ed2   },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL   },
};

int main (int argc, char *argv[])
{
  struct chr_block *chr_block_head1 = NULL; //for peak1
  struct chr_block *chr_block_head2 = NULL; //for peak2

  struct output *output_head = NULL; //for output
  struct output *ov_head = NULL; //for overlapping peaks
  struct output *nonov_head = NULL; //for non-overlapping peaks

  struct chr_block *ch1, *ch2; //for "for loop of chr"
  struct bs *bs_nonov; //for bs of file1 of which chr is not in file2

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
  if (file1 == NULL || file2 == NULL) usage();

  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
Input file1:                     %s\n\
Input file2:                     %s\n\
File1 column of chr, start, end: %d, %d, %d\n\
File2 column of chr, start, end: %d, %d, %d\n\
header flag:                     %d\n\
time:                            %s\n",\
 "ga_overlap", file1, file2, col_chr1, col_st1, col_ed1, col_chr2, col_st2, col_ed2, hf, ctime(&timer) );

  ga_parse_chr_bs(file1, &chr_block_head1, col_chr1, col_st1, col_ed1, hf); //parsing each binding sites for each chromosome
  ga_parse_chr_bs(file2, &chr_block_head2, col_chr2, col_st2, col_ed2, hf);

  for (ch1 = chr_block_head1; ch1; ch1 = ch1->next) {
    for (ch2 = chr_block_head2; ch2; ch2 = ch2->next) {
      if (!strcmp(ch1->chr, ch2->chr)) { //comparing peaks on the same chr
        cmp_overlap (ch1->bs_list, ch2->bs_list, &output_head, &ov_head, &nonov_head);
        break;
      }
    }
    if (ch2 == NULL) { //if chr_block2 doesn't have ch1
      for (bs_nonov = ch1->bs_list; bs_nonov; bs_nonov = bs_nonov->next) {
        ga_output_add (&nonov_head, bs_nonov->line);
        if (add_one_val(ga_line_out, bs_nonov->line, "NonOver\n") < 0) goto err; //making output link list with flag
        ga_output_add(&output_head, ga_line_out);
      }
    }
  }

  ga_parse_file_path (file1, path1, fn1, ext1); //parsing input file name into path, file name, and extension
  ga_parse_file_path (file2, path2, fn2, ext2);

  sprintf(output_name, "%s%s_over_%s%s", path1, fn1, fn2, ext1); //concatenating output file name
  ga_write_lines (output_name, ov_head, ga_header_line); //writing overlapping peaks

  sprintf(output_name, "%s%s_nonover_%s%s", path1, fn1, fn2, ext1); //concatenating output file name
  ga_write_lines (output_name, nonov_head, ga_header_line); //writing non-overlapping peaks

  sprintf(output_name, "%s%s_vs_%s%s", path1, fn1, fn2, ext1);
  if (ga_header_line != NULL) { //if header line
    if (add_one_val(ga_line_out, ga_header_line, "overlap_flag\n") < 0) goto err; //adding one extra column
    ga_write_lines (output_name, output_head, ga_line_out); //note that header is line_out, not ga_header_line
  }
  else ga_write_lines (output_name, output_head, ga_header_line);

  if (ga_header_line) free(ga_header_line);
  ga_free_chr_block(&chr_block_head1);
  ga_free_chr_block(&chr_block_head2);

  ga_free_output(&output_head);
  ga_free_output(&ov_head);
  ga_free_output(&nonov_head);

  return 0;

err:
  return -1;
}

/*
 * This function compares overlapping of bs.
 * *bs1: pointer to bs link list for file1
 * *bs2: pointer to bs link list for file2
 * **output_head: pointer of pointer to output link list for output file with flag
 * **ov_head: pointer of pointer to output link list for output file overlapping bs2
 * **nonov_head: pointer of pointer to output link list for output file not overlapping bs2
 */
static void cmp_overlap (struct bs *bs1, struct bs *bs2, struct output **output_head, struct output **ov_head, struct output **nonov_head)
{
//  char line_out[LINE_STR_LEN];
  struct bs *i, *j;

  for (i = bs1; i; i = i->next) {
    for (j = bs2; j; j = j->next) {
      if (i->ed >= j->st && i->st <= j->ed) { //checking the overlapping
        ga_output_add(ov_head, i->line);
        if (add_one_val(ga_line_out, i->line, "Over\n") < 0) goto err;
        ga_output_add(output_head, ga_line_out);
        break;
      }
    }
    if (j == NULL) { //if i is not overlapped with any bs2
      ga_output_add(nonov_head, i->line);
      if (add_one_val(ga_line_out, i->line, "NonOver\n") < 0) goto err;
      ga_output_add(output_head, ga_line_out);
    }
  }

err:
  return;
}

/*
 * This function add one more value to string. If line_out[xxx], line = "aaa\tbbb\n", val = "ccc\n", line_out is "aaa\tbbb\tccc\n".
 * line_out[]: char array. This must have size of LINE_STR_LEN.
 * *line: pointer to char to be added.
 * *val: pointer to char for adding. Put '\n' at the last position if you need.
 */
static int add_one_val (char line_out[], const char *line, const char *val)
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
