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

static int cmp_overlap (struct bs *bs1, struct bs *bs2, struct output **output_head, struct output **ov_head, struct output **nonov_head);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_overlap\n\
Summary: report the overlaps of two peaks/summits\n\n\
Usage:   ga_overlap [options] -1 <file1> -2 <file2>\n\
   or:   ga_overlap [options] -1 <file1> -2 <file2> --hw <half_window: int> --col_start1 <int> --col_end1 <int> ##col_start1 and col_end1 should be same for -hw option\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the header of file1 is preserved.\n\
         --count: report the number of overlapped peaks of file2 for each peak of file1.\n\
         --hw <int>: the peak width is extended for this fixed half window from summit of file1.\n\
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
  printf("'ga_overlap' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}

static int hf = 0;
static int cf = 0;
static int hw = 0;
static char *file1 = NULL;
static char *file2 = NULL;
static int col_chr1 = 0, col_chr2 = 0;
static int col_st1 = 1, col_st2 = 1;
static int col_ed1 = 2, col_ed2 = 2;
char *ga_header_line = NULL; //header line. Note this is external global variable
static char ga_line_out[LINE_STR_LEN]; //output line with overlapping flag
static double totnb = 0.0, ovnb = 0.0, novnb = 0.0; //peak numbers


static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage    },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage    },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version  },
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf      },
  {"--count"     , ARGUMENT_TYPE_FLAG_ON , &cf      },
  {"--hw"        , ARGUMENT_TYPE_INTEGER , &hw      },
  {"-1"          , ARGUMENT_TYPE_STRING  , &file1   },
  {"-2"          , ARGUMENT_TYPE_STRING  , &file2   },
  {"--col_chr1"  , ARGUMENT_TYPE_INTEGER , &col_chr1},
  {"--col_chr2"  , ARGUMENT_TYPE_INTEGER , &col_chr2},
  {"--col_start1", ARGUMENT_TYPE_INTEGER , &col_st1 },
  {"--col_start2", ARGUMENT_TYPE_INTEGER , &col_st2 },
  {"--col_end1"  , ARGUMENT_TYPE_INTEGER , &col_ed1 },
  {"--col_end2"  , ARGUMENT_TYPE_INTEGER , &col_ed2 },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL     },
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

  FILE *fp;

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
counter flag:                    %d\n\
half window:                     %d\n\
time:                            %s\n",\
 "ga_overlap", file1, file2, col_chr1, col_st1, col_ed1, col_chr2, col_st2, col_ed2, hf, cf, hw, ctime(&timer) );

  ga_parse_chr_bs(file1, &chr_block_head1, col_chr1, col_st1, col_ed1, -1, hf); //parsing each binding sites for each chromosome without strand info
  ga_parse_chr_bs(file2, &chr_block_head2, col_chr2, col_st2, col_ed2, -1, hf);

  for (ch1 = chr_block_head1; ch1; ch1 = ch1->next) {
    for (ch2 = chr_block_head2; ch2; ch2 = ch2->next) {
      if (!strcmp(ch1->chr, ch2->chr)) { //comparing peaks on the same chr
        if(cmp_overlap (ch1->bs_list, ch2->bs_list, &output_head, &ov_head, &nonov_head) != 0){
          LOG("error: error in cmp_overlap function.");
          goto err;
        }
        break;
      }
    }
    if (ch2 == NULL) { //if chr_block2 doesn't have ch1
      for (bs_nonov = ch1->bs_list; bs_nonov; bs_nonov = bs_nonov->next) {
        ga_output_add (&nonov_head, bs_nonov->line); //caution: the order is reversed 
        if (cf) { //if count
          if (add_one_val(ga_line_out, bs_nonov->line, "NonOver\t0\n") < 0) {
            LOG("error: output line was too long.");
            goto err; //making output link list with flag
          }
        } else {
          if (add_one_val(ga_line_out, bs_nonov->line, "NonOver\n") < 0) {
            LOG("error: output line was too long.");
            goto err; //making output link list with flag
          }
        }
        ga_output_add(&output_head, ga_line_out); //caution: the order is reversed
        totnb += 1.0; //total peak number
        novnb += 1.0; //non-overlapping peak number
      }
    }
  }

  ga_parse_file_path (file1, path1, fn1, ext1); //parsing input file name into path, file name, and extension
  ga_parse_file_path (file2, path2, fn2, ext2);

  if (hw) sprintf(output_name, "%s%s_hw%d_over_%s%s", path1, fn1, hw, fn2, ext1); //concatenating output file name
  else sprintf(output_name, "%s%s_over_%s%s", path1, fn1, fn2, ext1); //concatenating output file name
  ga_write_lines (output_name, ov_head, ga_header_line); //writing overlapping peaks

  if (hw) sprintf(output_name, "%s%s_hw%d_nonover_%s%s", path1, fn1, hw, fn2, ext1); //concatenating output file name
  else sprintf(output_name, "%s%s_nonover_%s%s", path1, fn1, fn2, ext1); //concatenating output file name
  ga_write_lines (output_name, nonov_head, ga_header_line); //writing non-overlapping peaks

  if (hw) sprintf(output_name, "%s%s_hw%d_vs_%s%s", path1, fn1, hw, fn2, ext1);
  else sprintf(output_name, "%s%s_vs_%s%s", path1, fn1, fn2, ext1);
  if (ga_header_line != NULL) { //if header line
    if (cf) {
      if (add_one_val(ga_line_out, ga_header_line, "overlap_flag\tcount\n") < 0) {
        LOG("error: output line was too long.");
        goto err; //adding one extra column
      }
    } else {
      if (add_one_val(ga_line_out, ga_header_line, "overlap_flag\n") < 0) {
        LOG("error: output line was too long.");
        goto err; //adding one extra column
      }
    }
    ga_write_lines (output_name, output_head, ga_line_out); //note that header is line_out, not ga_header_line
  }
  else ga_write_lines (output_name, output_head, ga_header_line);

  if (hw) sprintf(output_name, "%s%s_hw%d_vs_%s_summary.txt", path1, fn1, hw, fn2); //concatenating output file name
  else sprintf(output_name, "%s%s_vs_%s_summary.txt", path1, fn1, fn2); //concatenating output file name
  if ((fp = fopen(output_name, "w")) == NULL) {
    LOG("error: output file cannot be open.");
    goto err;
  }
  sprintf(ga_line_out, "name\ttotal\tOver\tNonOver\n");
  if (fputs(ga_line_out, fp) == EOF) {
    LOG("error: file writing error.");
    goto err;
  }
  sprintf(ga_line_out, "%s\t%d\t%d(%.3f)\t%d(%.3f)\n", fn1, (int)totnb, (int)ovnb, ovnb/totnb, (int)novnb, novnb/totnb);
  if (fputs(ga_line_out, fp) == EOF) {
    LOG("error: file writing error.");
    goto err;
  }
  fclose(fp);

  if (ga_header_line) free(ga_header_line);
  ga_free_chr_block(&chr_block_head1);
  ga_free_chr_block(&chr_block_head2);

  ga_free_output(&output_head);
  ga_free_output(&ov_head);
  ga_free_output(&nonov_head);

  return 0;

err:
  if (ga_header_line) free(ga_header_line);
  if (chr_block_head1) ga_free_chr_block(&chr_block_head1);
  if (chr_block_head2) ga_free_chr_block(&chr_block_head2);

  if (output_head) ga_free_output(&output_head);
  if (ov_head) ga_free_output(&ov_head);
  if (nonov_head) ga_free_output(&nonov_head);
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
static int cmp_overlap (struct bs *bs1, struct bs *bs2, struct output **output_head, struct output **ov_head, struct output **nonov_head)
{
  struct bs *i, *j;
  char tmp[64];
  int st1, ed1; //start and end position
  int c; //counter

  for (i = bs1; i; i = i->next) {
    c = 0; //init counter
    if (hw) {
      st1 = i->st - hw; //extending the peak width
      ed1 = i->st + hw; //extending the peak width
    } else {
      st1 = i->st;
      ed1 = i->ed;
    }

    for (j = bs2; j; j = j->next) {
//      if (i->ed >= j->st && i->st <= j->ed) { //checking the overlapping
      if (ed1 >= j->st && st1 <= j->ed) { //checking the overlapping
        if (!c) { //ov_head is added only once
          ga_output_add(ov_head, i->line); //caution: the order is reversed
          ovnb += 1.0; //overlapping peak number
        }
        c++;
        if(!cf) {
          if (add_one_val(ga_line_out, i->line, "Over\n") < 0) {
            LOG("error: output line is too long.");
            goto err;
          }
          ga_output_add(output_head, ga_line_out); //caution: the order is reversed
          break;
        }
      }
    }
    if (cf && c) { //if at least one peak was overlapped for i
      sprintf(tmp, "Over\t%d\n", c);
      if (add_one_val(ga_line_out, i->line, tmp) < 0) {
        LOG("error: output line is too long.");
        goto err;
      }
      ga_output_add(output_head, ga_line_out); //caution: the order is reversed
    }
//    if (j == NULL) { //if i is not overlapped with any bs2
    if (!c) { //if i is not overlapped with any bs2
      ga_output_add(nonov_head, i->line); //caution: the order is reversed
      novnb += 1.0; //non-overlapping peak number
      if (cf) {
        if (add_one_val(ga_line_out, i->line, "NonOver\t0\n") < 0) {
          LOG("error: output line is too long.");
          goto err;
        }
        ga_output_add(output_head, ga_line_out); //caution: the order is reversed
      } else {
        if (add_one_val(ga_line_out, i->line, "NonOver\n") < 0) {
          LOG("error: output line is too long.");
          goto err;
        }
        ga_output_add(output_head, ga_line_out); //caution: the order is reversed
      }
    }
    totnb += 1.0; //total peak number
  }
  return 0;

err:
  return -1;
}


