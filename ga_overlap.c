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
#include "ga_my.h"

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
  printf("Tool:    ga_overlap\n\n\
Summary: report the overlaps of two peaks/summits\n\n\
Usage:   ga_overlap [options] -1 <file1> -2 <file2>\n\
   or:   ga_overlap [options] -1 <file1> -2 <file2> --hw <half_window: int> --col_start1 <int> --col_end1 <int> ##col_start1 and col_end1 should be same for -hw option\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the header of file1 is preserved (default: off).\n\
         --count: report the number of overlapped peaks of file2 for each peak of file1 (default: off).\n\
         --preserve2: preserve file2 content if overlapped (default: off).\n\
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
static char hfs[4] = "off\0";
static int cf = 0;
static char cfs[4] = "off\0";
static int p2 = 0;
static char p2s[4] = "off\0";
static int hw = 0;
static char *file1 = NULL;
static char *file2 = NULL;
static int col_chr1 = 0, col_chr2 = 0;
static int col_st1 = 1, col_st2 = 1;
static int col_ed1 = 2, col_ed2 = 2;
char *ga_header_line = NULL; //header line. Note this is external global variable
static char ga_line_out[LINE_STR_LEN] = {0}; //output line with overlapping flag
static char line2[LINE_STR_LEN] = {0}; //line for file2
static double totnb = 0.0, ovnb = 0.0, novnb = 0.0; //peak numbers


static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage    },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage    },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version  },
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf      },
  {"--count"     , ARGUMENT_TYPE_FLAG_ON , &cf      },
  {"--preserve2" , ARGUMENT_TYPE_FLAG_ON , &p2      },
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
  char path1[PATH_STR_LEN] = {0};
  char fn1[FILE_STR_LEN] = {0};
  char ext1[EXT_STR_LEN] = {0};
  char path2[PATH_STR_LEN] = {0};
  char fn2[FILE_STR_LEN] = {0};
  char ext2[EXT_STR_LEN] = {0};
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN] = {0}; //output file name
  char *ga_header_line2 = NULL; //header line2 if p2
  int i=0;
  time_t timer;

  FILE *fp;

  argument_read(&argc, argv, args);//reading arguments
  if (file1 == NULL || file2 == NULL) usage();

  if(hf) strcpy(hfs, "on\0");
  if(cf) strcpy(cfs, "on\0");
  if(p2) strcpy(p2s, "on\0");
  time(&timer);
  printf("Tool:                            %s\n\n\
Input file1:                     %s\n\
Input file2:                     %s\n\
File1 column of chr, start, end: %d, %d, %d\n\
File2 column of chr, start, end: %d, %d, %d\n\
header flag:                     %s\n\
counter flag:                    %s\n\
preserve file2?:                 %s\n\
half window:                     %d\n\
time:                            %s\n",\
 "ga_overlap", file1, file2, col_chr1, col_st1, col_ed1, col_chr2, col_st2, col_ed2, hfs, cfs, p2s, hw, ctime(&timer) );

  ga_parse_chr_bs(file2, &chr_block_head2, col_chr2, col_st2, col_ed2, -1, hf);
  if (p2 && NULL != ga_header_line) ga_header_line2 = strdup(ga_header_line); //preserving header2 
  if (NULL != ga_header_line) my_free(ga_header_line);
  ga_header_line = NULL;
  ga_parse_chr_bs(file1, &chr_block_head1, col_chr1, col_st1, col_ed1, -1, hf); //parsing each binding sites for each chromosome without strand info

  //making "tab" line for non-overlapping...
  if (p2) {
    bs_nonov = chr_block_head2 -> bs_list; //the very first line
    while(bs_nonov -> line [i] != '\0'){
      if (bs_nonov -> line [i] == '\t'){
        if (strlen(line2) + strlen("NA\t") + 1 < LINE_STR_LEN) strncat(line2, "NA\t", strlen("NA\t"));
      }
      i++;
    }
    if (strlen(line2) + strlen("NA\t") + 1 < LINE_STR_LEN) strncat(line2, "NA\t", strlen("NA\t"));
  }
  //

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

        sprintf(ga_line_out, "%s", bs_nonov->line);
        ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t

        if (p2) {
          if (strlen(ga_line_out) + strlen(line2) + 1 < LINE_STR_LEN) strncat(ga_line_out, line2, strlen(line2));
        }

        if (strlen(ga_line_out) + strlen("NonOver\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "NonOver\t", strlen("NonOver\t"));

        if (cf) {
          if (strlen(ga_line_out) + strlen("0\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "0\t", strlen("0\t"));
        }

        ga_line_out[strlen(ga_line_out) - 1] = '\n'; //\t was replaced by \n
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
    sprintf(ga_line_out, "%s", ga_header_line);
    ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t

    if (p2) {
      if (NULL == ga_header_line2) {
        LOG("error: only file1 had a header.");
        goto err;
      }
      if (strlen(ga_line_out) + strlen(ga_header_line2) + 1 < LINE_STR_LEN) strncat(ga_line_out, ga_header_line2, strlen(ga_header_line2));
      ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t
    }

    if (strlen(ga_line_out) + strlen("overlap_flag\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "overlap_flag\t", strlen("overlap_flag\t"));

    if (cf) {
      if (strlen(ga_line_out) + strlen("count\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "count\t", strlen("count\t"));
    }
    ga_line_out[strlen(ga_line_out) - 1] = '\n'; //\t was replaced by \n
    ga_write_lines (output_name, output_head, ga_line_out);
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

  if (ga_header_line) my_free(ga_header_line);
  if (ga_header_line2) my_free(ga_header_line2);
  ga_free_chr_block(&chr_block_head1);
  ga_free_chr_block(&chr_block_head2);

  ga_free_output(&output_head);
  ga_free_output(&ov_head);
  ga_free_output(&nonov_head);

  return 0;

err:
  if (ga_header_line) my_free(ga_header_line);
  if (ga_header_line2) my_free(ga_header_line2);
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
  struct bs *i, *j, *j_tmp;
  char tmp[80] = {0};
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
      if (ed1 >= j->st && st1 <= j->ed) { //checking the overlapping
        if (!c) { //ov_head is added only once
          ga_output_add(ov_head, i->line); //caution: the order is reversed
          ovnb += 1.0; //overlapping peak number
        }
        c++;
        if(!cf) {
          sprintf(ga_line_out, "%s", i->line);
          ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t

          if (p2) {
            if (strlen(ga_line_out) + strlen(j->line) + 1 < LINE_STR_LEN) strncat(ga_line_out, j->line, strlen(j->line));
            ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t
          }

          if (strlen(ga_line_out) + strlen("Over\n") + 1 < LINE_STR_LEN) strncat(ga_line_out, "Over\n", strlen("Over\n"));

          ga_output_add(output_head, ga_line_out); //caution: the order is reversed
          break;
        } //if !cf
        j_tmp = j; //preserving overlapping line2
      } //if overlap
    } //for j
    if (cf && c) { //if at least one peak was overlapped for i
      sprintf(ga_line_out, "%s", i->line);
      ga_line_out[strlen(ga_line_out) - 1] = '\t';

      if (p2) {
        if (strlen(ga_line_out) + strlen(j_tmp->line) + 1 < LINE_STR_LEN) strncat(ga_line_out, j_tmp->line, strlen(j_tmp->line));
        ga_line_out[strlen(ga_line_out) - 1] = '\t';
      }

      sprintf(tmp, "Over\t%d\n", c);
      if (strlen(ga_line_out) + strlen(tmp) + 1 < LINE_STR_LEN) strncat(ga_line_out, tmp, strlen(tmp));

      ga_output_add(output_head, ga_line_out); //caution: the order is reversed
    } //if cf && c
    if (!c) { //if i is not overlapped with any bs2
      ga_output_add(nonov_head, i->line); //caution: the order is reversed
      novnb += 1.0; //non-overlapping peak number

      sprintf(ga_line_out, "%s", i->line);
      ga_line_out[strlen(ga_line_out) - 1] = '\t'; //\n was replaced by \t

      if (p2) {
        if (strlen(ga_line_out) + strlen(line2) + 1 < LINE_STR_LEN) strncat(ga_line_out, line2, strlen(line2));
      }

      if (strlen(ga_line_out) + strlen("NonOver\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "NonOver\t", strlen("NonOver\t"));

      if (cf) {
        if (strlen(ga_line_out) + strlen("0\t") + 1 < LINE_STR_LEN) strncat(ga_line_out, "0\t", strlen("0\t"));
      }
      ga_line_out[strlen(ga_line_out) - 1] = '\n'; //\t was replaced by \n
      ga_output_add(output_head, ga_line_out); //caution: the order is reversed
    }
    totnb += 1.0; //total peak number
  }
  return 0;
}


