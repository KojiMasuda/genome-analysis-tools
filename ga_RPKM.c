/*
 * This program is one of the genome analysis tools.
 * This program calculates expression level as reads per kilobase per million mapped reads (RPKM). Input file is normalized mapped reads (per million) in bedgraph format.
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "sort_list.h"
#include "argument.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)


static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_RPKM\n\
Summary: report the expression level as reads per kilobase per million mapped reads (RPKM)\n\n\
Usage:   ga_RPKM [options] -exp <expression file in bedgraph format> -ref <reference file> -readlen <read length(default:101)>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --col_chr <int>: column number for chromosome of reference file (default:2).\n\
         --col_start <int>: column number for gene start position of reference file (default:4).\n\
         --col_end <int>: column number for gene end position of reference file (default:5).\n\
         --col_strand <int>: column number for strand of reference file (default:3).\n\
         --col_exon_st <int>: column number for exon start positions of reference file (default:9).\n\
         --col_exon_ed <int>: column number for exon end positions of reference file (default:10).\n\
         --header: the first line of reference file is considered as header (default:off).\n");
  exit(0);
}

static void version()
{
  printf("'ga_RPKM' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}


static int hf = 0;
static char *fileexp = NULL;
static char *fileref = NULL;
static int readlen = 101;
static int col_chr = 2;
static int col_st = 4;
static int col_ed = 5;
static int col_strand = 3;
static int col_ex_st = 9;
static int col_ex_ed = 10;
char *ga_header_line = NULL; //header line. Note this is external global variable
static char ga_line_out[LINE_STR_LEN] = {0}; //output line

static const Argument args[] = {
  {"-h"           , ARGUMENT_TYPE_FUNCTION, usage        },
  {"--help"       , ARGUMENT_TYPE_FUNCTION, usage        },
  {"-v"           , ARGUMENT_TYPE_FUNCTION, version      },
  {"--header"     , ARGUMENT_TYPE_FLAG_ON , &hf          },
  {"-exp"         , ARGUMENT_TYPE_STRING  , &fileexp     },
  {"-ref"         , ARGUMENT_TYPE_STRING  , &fileref     },
  {"-readlen"     , ARGUMENT_TYPE_INTEGER , &readlen     },
  {"--col_chr"    , ARGUMENT_TYPE_INTEGER , &col_chr     },
  {"--col_start"  , ARGUMENT_TYPE_INTEGER , &col_st      },
  {"--col_end"    , ARGUMENT_TYPE_INTEGER , &col_ed      },
  {"--col_strand" , ARGUMENT_TYPE_INTEGER , &col_strand  },
  {"--col_exon_st", ARGUMENT_TYPE_INTEGER , &col_ex_st   },
  {"--col_exon_ed", ARGUMENT_TYPE_INTEGER , &col_ex_ed   },
  {NULL           , ARGUMENT_TYPE_NONE    , NULL         },
};


int main (int argc, char *argv[])
{
  argument_read(&argc, argv, args);//reading arguments
  if (fileexp == NULL || fileref == NULL) usage();

  struct chr_block *chr_block_headexp = NULL; //for expression
  struct chr_block *chr_block_headref = NULL; //for reference

  struct chr_block *ch_ref, *ch_exp; //for "for loop of chr"
  struct ref *r;
  struct sig *j1, *j1_tmp;
  struct output *output_head = NULL; //for output

  /*path, filename, and extension*/
  char path_ref[PATH_STR_LEN] = {0};
  char fn_ref[FILE_STR_LEN] = {0};
  char ext_ref[EXT_STR_LEN] = {0};
  char path_exp[PATH_STR_LEN] = {0};
  char fn_exp[FILE_STR_LEN] = {0};
  char ext_exp[EXT_STR_LEN] = {0};
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN] = {0}; //output file name
  char tmp[LINE_STR_LEN] = {0};
  char exon_st_tmp[10000] = {0}, exon_ed_tmp[10000] = {0}; //exon positions, like "123,456,"
  char exon_st[80] = {0}, exon_ed[80] = {0}; //exon position like "123"
  char *end; //for strtoul
  unsigned long arr_ex_st[1000]; //array of exon positions
  unsigned long arr_ex_ed[1000];
  unsigned long st;
  unsigned long tmp_st;
  unsigned long ed;
  unsigned long tmp_ed;
  unsigned long exon_len;
  int fl, c, e, e_tmp, exon_nb;
  double val;

  time_t timer;


  time(&timer);
  printf("Program:                           %s\n\
Tools:                             genome analysis tools\n\n\
Input expression file:             %s\n\
Input reference file:              %s\n\
Read length:                       %d\n\
ref col of chr, start, end, strand:%d, %d, %d, %d\n\
ref col exon start, end:           %d, %d\n\
header flag:                       %d\n\
time:                              %s\n",\
 "ga_RPKM", fileexp, fileref, readlen, col_chr, col_st, col_ed, col_strand, col_ex_st, col_ex_ed, hf, ctime(&timer) );

  ga_parse_file_path (fileexp, path_exp, fn_exp, ext_exp); //parsing input file name into path, file name, and extension
  ga_parse_file_path (fileref, path_ref, fn_ref, ext_ref);

  ga_parse_chr_ref (fileref, &chr_block_headref, col_chr, col_st, col_ed, col_strand, col_ex_st, col_ex_ed, hf); //parsing each ref for each chromosome

  // reading exp file
  ga_parse_bedgraph (fileexp, &chr_block_headexp);

  // sorting summit and sig
  chr_block_headref = ga_mergesort_chr(chr_block_headref);
  chr_block_headexp = ga_mergesort_chr(chr_block_headexp);

  for (ch_ref = chr_block_headref; ch_ref; ch_ref = ch_ref -> next) {
    ch_ref -> ref_list = ga_mergesort_ref(ch_ref -> ref_list);
  }
  for (ch_exp = chr_block_headexp; ch_exp; ch_exp = ch_exp -> next) {
    ch_exp -> sig_list = ga_mergesort_sig(ch_exp -> sig_list);
  }

  for (ch_ref = chr_block_headref; ch_ref; ch_ref = ch_ref -> next) {
    for (ch_exp = chr_block_headexp; ch_exp; ch_exp = ch_exp -> next) {
      if (!strcmp(ch_ref->chr, ch_exp->chr)) break;
    }

    if (ch_exp == NULL) {
      for (r = ch_ref -> ref_list; r; r = r -> next) {
        sprintf(tmp, "0\n");
        if (add_one_val(ga_line_out, r -> line, tmp) < 0) {
          LOG("error: output line was too long.");
          goto err;
        }
        ga_output_append (&output_head, ga_line_out);
      }
    }

    j1_tmp = NULL; //the "marker" of signal position to speed up the calc. j1_tmp is the left most position for each ref.
    for (r = ch_ref -> ref_list; r; r = r -> next) {
      fl = 0; //init of flag
      e = 0; //char position for exon
      e_tmp = 0; //char position for tmp exon
      c = 0; //number of exon
      memset(exon_st_tmp, '\0', sizeof(exon_st_tmp));
      memset(exon_ed_tmp, '\0', sizeof(exon_ed_tmp));
      strncpy(exon_st_tmp, r -> ex_st, sizeof(exon_st_tmp) - 1); //exon st positions, like "123,456,"
      strncpy(exon_ed_tmp, r -> ex_ed, sizeof(exon_ed_tmp) - 1); //exon ed positions, like "234,567,"
      val = 0;
      exon_len = 0;

      while(exon_st_tmp[e_tmp] != '\0') {
        if (exon_st_tmp[e_tmp] == ',') {
          exon_st[e] = '\0';
          arr_ex_st[c] = strtoul(exon_st, &end, 0);
          e = 0;
          e_tmp++;
          c++;
          continue;
        } else {
          exon_st[e] = exon_st_tmp[e_tmp]; //extracting start pos of each exon
        }
        e++;
        e_tmp++;
      } //exon start

      e = 0;
      e_tmp = 0;
      c = 0;
      while(exon_ed_tmp[e_tmp] != '\0') {
        if (exon_ed_tmp[e_tmp] == ',') {
          exon_ed[e] = '\0';
          arr_ex_ed[c] = strtoul(exon_ed, &end, 0);
          e = 0;
          e_tmp++;
          c++;
          continue;
        } else {
          exon_ed[e] = exon_ed_tmp[e_tmp]; //extracting end pos of each exon
        }
        e++;
        e_tmp++;
      } //exon end
      exon_nb = c; //exon number

      for (c = 0; c < exon_nb; c++) {
        st = arr_ex_st[c] - 1; //-1 because bedgraph is zero-based, exon st also should be zero-based
        ed = arr_ex_ed[c]; //but end is half-open, so this is not -1.
        exon_len += ed - st; 
        if (j1_tmp == NULL) {
          for (j1 = ch_exp->sig_list; j1; j1=j1->next) {
            if (st < j1->ed && j1->st < ed) {
              break; //if one of sig block is inside the exon
            } else if (j1->st >= ed) { //if there's no chance for j1 to overlap exon
              j1 = NULL;
              break;
            }
          }
        } else {
          for (j1 = j1_tmp; j1; j1=j1->next) { //starting the search from j1_tmp to speed up!
            if (st < j1->ed && j1->st < ed) {
              break; //if one of sig block is inside the exon
            } else if (j1->st >= ed) { //if there's no chance for j1 to overlap exon
              j1 = NULL;
              break;
            }
          }
        }

        if (j1 == NULL) continue; //if the exon doesn't overlap with any sig

        if (j1_tmp == NULL) { //if j1_tmp is not set for the chr
          j1_tmp = j1;
          fl = 1;
        } else if (!fl) { //else if j1_tmp is not set for the ref (fl == 0)
          j1_tmp = j1;
          fl = 1;
        }

        for (;j1 ; j1 = j1->next) {
          if (j1->st >= ed) break; //if the sig pos is out of the win
          if (st > j1->st) tmp_st = st; //if st of sig block is left-side to pos of st
          else tmp_st = j1->st;
          if (j1->ed > ed) tmp_ed = ed; //if ed of sig block is right-side to pos of ed
          else tmp_ed = j1->ed;
          val += (j1->val) * (tmp_ed - tmp_st);
        }
      } //for exon

      sprintf(tmp, "%.6f\n", val / ( (float)exon_len / 1000 ) / (float)readlen ); //calculating RPKM
      if (add_one_val(ga_line_out, r -> line, tmp) < 0) {
        LOG("error: output line was too long.");
        goto err;
      }
      ga_output_append (&output_head, ga_line_out);
    } //for ref
  } //for chr

  sprintf(output_name, "%sRPKM_of_%s_%s.txt", path_exp, fn_exp, fn_ref);
  ga_write_lines (output_name, output_head, ga_header_line);

  if (chr_block_headref) ga_free_chr_block(&chr_block_headref);
  if (chr_block_headexp) ga_free_chr_block(&chr_block_headexp);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);
  return 0;

err:
  if (chr_block_headref) ga_free_chr_block(&chr_block_headref);
  if (chr_block_headexp) ga_free_chr_block(&chr_block_headexp);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);
  return -1;
}
