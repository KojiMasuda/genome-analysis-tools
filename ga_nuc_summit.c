/*
 * This program is one of the genome analysis tools.
 * This program calculte nucleotide content distributions around summits.
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


static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_nuc_summit\n\
Summary: report the nucleotide composition(number of ATCG) or AT content around summit.\n\n\
Usage:   ga_nuc_summit -fa <fasta> -smt <summit file> -gt <genome table> -n_flag <nucleotide flag:A | T | C | G | AT>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --header: the first line of the summit file is considered as header (default:off).\n\
         --col_chr: <int> column number for chromosome of summit file (default:0).\n\
         --col_start: <int> column number for start position of summit file (default:1).\n\
         --col_end: <int> column number for end position of summit file (default:2).\n\
         --col_strand: <int> column number for strand of summit file (default:-1).\n\
         --hw: <int> half range size (default:1000).\n\
         --step: <int> step size (default:10).\n\
         --win: <int> window size (default:100).\n");
  exit(0);
}

static void version()
{
  printf("'ga_nuc_summit' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}

char *ga_header_line = NULL; //header line. Note this is external global variable
static char *smt = NULL; //summit file
static char *fa = NULL; //fasta
static char *gt = NULL; //genome table
static char *n_flag = NULL; //nuc flag
static char ga_line_out[LINE_STR_LEN] = {0}; //output line with nucleotide composition
static int col_chr = 0;
static int col_st = 1;
static int col_ed = 2;
static int col_strand = -1;
static int hf = 0;
static int hw = 1000;
static int step = 10;
static int win = 100;

static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage       },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage       },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version     },
  {"-smt"        , ARGUMENT_TYPE_STRING  , &smt        },
  {"-fa"         , ARGUMENT_TYPE_STRING  , &fa         },
  {"-gt"         , ARGUMENT_TYPE_STRING  , &gt         },
  {"-n_flag"     , ARGUMENT_TYPE_STRING  , &n_flag     },
  {"--header"    , ARGUMENT_TYPE_FLAG_ON , &hf         },
  {"--col_chr"   , ARGUMENT_TYPE_INTEGER , &col_chr    },
  {"--col_start" , ARGUMENT_TYPE_INTEGER , &col_st     },
  {"--col_end"   , ARGUMENT_TYPE_INTEGER , &col_ed     },
  {"--col_strand", ARGUMENT_TYPE_INTEGER , &col_strand },
  {"--hw"        , ARGUMENT_TYPE_INTEGER , &hw         },
  {"--step"      , ARGUMENT_TYPE_INTEGER , &step       },
  {"--win"       , ARGUMENT_TYPE_INTEGER , &win        },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL        },
};

int main (int argc, char *argv[])
{
  /*path, filename, and extension*/
  char path[PATH_STR_LEN] = {0};
  char fn[FILE_STR_LEN] = {0};
  char ext[EXT_STR_LEN] = {0};
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN] = {0}; //output file name
  char tmp[80] = {0}; //the value for each window
  char *frag; //fragment

  int winNb, w; //window number

  struct chr_block *chr_block_head_smt = NULL; //for parsing summit file
  struct chr_block *ch1 ; //for 
  struct chr_block *chr_block_head_gt = NULL; //for parsing genome table
  struct chr_block_fa *chr_block_head_fa = NULL; //for parsing fasta
  struct chr_block_fa *ch2; //for
  struct bs *bs; //for
  struct output *out_head; //output header

  unsigned long i, c_A, c_T, c_G, c_C; //count of ATCG
  long s; //start position of the window (can be negative if summit is close to position 0 and hw is too large...)

  time_t timer;

  argument_read(&argc, argv, args);//reading arguments
  if (fa == NULL || smt == NULL || gt == NULL) usage();

  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
fasta file:                      %s\n\
summit file:                     %s\n\
genome table:                    %s\n\
nucleotide flag:                 %s\n\
header in summit file?:          %d\n\
col for chr, start, end:         %d, %d, %d\n\
col for strand?:                 %d\n\
half range:                      %d\n\
step size:                       %d\n\
window size:                     %d\n\
time:                            %s\n",\
 "ga_nuc_smt", fa, smt, gt, n_flag, hf, col_chr, col_st, col_ed, col_strand, hw, step, win, ctime(&timer) );

  ga_parse_chr_bs (smt, &chr_block_head_smt, col_chr, col_st, col_ed, col_strand, hf); //parsing summit file
  chr_block_head_smt = ga_mergesort_chr(chr_block_head_smt); // sorting summit chr

  for (ch1 = chr_block_head_smt; ch1; ch1 = ch1 -> next) {
    ch1 -> bs_list = ga_mergesort_bs(ch1 -> bs_list); //sorting bs
  }

  ga_parse_chr_bs (gt, &chr_block_head_gt, 0, 1, 1, -1, 0); //parsing genome tible
  if (ga_parse_chr_fa(fa, &chr_block_head_fa, chr_block_head_gt) != 0){ //parsing fasta
    LOG("error: error in ga_parse_chr_fa function.");
    goto err;
  }
  ga_parse_file_path (smt, path, fn, ext); //parsing input file name into path, file name, and extension

  frag = (char *)malloc(sizeof(char) * (win + 1)); //allocating memory for fragment DNA from genome
  winNb = (2 * hw) / step + 1; //window number

  for (ch1 = chr_block_head_smt; ch1; ch1 = ch1->next) {
    for (ch2 = chr_block_head_fa; ch2; ch2 = ch2->next) {
      if (!(strcmp(ch1->chr, ch2->chr))) break;
    }

    if (ch2 == NULL) {
      printf("warning: chromosome %s was not found in the fasta file, skipped\n", ch1->chr);
      continue;
    }
    else printf("calculating on %s \n", ch1->chr);

    for (bs = ch1 -> bs_list; bs; bs = bs -> next) { //counting nucleotide for each summit
      sprintf(ga_line_out, "%s_%lu-%lu\t", ch1->chr, bs->st, bs->ed);
      if (bs->strand == '-') s = bs->ed + hw - win / 2 - 1; //if - strand
      else s = bs->st - hw - win / 2 - 1; //if + strand or no information

      for (w = 0; w < winNb; w++) {
        if (s < 0) { //if win position is less than 0.
          printf("warning: window %d of summit %lu on %s is less than zero. NA is created.\n", w+1, bs->st, ch1->chr);
          if (w == winNb - 1) sprintf(tmp, "NA\n");
          else sprintf(tmp, "NA\t");
          if (strlen(ga_line_out) + strlen(tmp) + 1 > sizeof(ga_line_out)) {
            LOG("error: string per line is too long.");
            goto err;
          }
          strcat(ga_line_out, tmp); //concatenating val for each window

          if (bs->strand == '-') s = s - step; //if - strand
          else s = s + step; //if + strand or no information
          continue;
        } else if (s + win > ch2->letter_len) { //if win position is over the letter
          printf("warning: window %d of summit %lu on %s is over the chromosome. NA is created.\n", w+1, bs->st, ch1->chr);
          if (w == winNb - 1) sprintf(tmp, "NA\n");
          else sprintf(tmp, "NA\t");
          if (strlen(ga_line_out) + strlen(tmp) + 1 > sizeof(ga_line_out)) {
            LOG("error: string per line is too long.");
            goto err;
          }
          strcat(ga_line_out, tmp); //concatenating val for each window

          if (bs->strand == '-') s = s - step; //if - strand
          else s = s + step; //if + strand or no information
          continue;
        }

        strncpy(frag, ch2->letter + s, win);  //sliced fragment DNA from genome. Starting from letter[0], win nucleotide is copied to frag. Next the slice starts from letter[0 + s].
        frag[win] = '\0'; //null
        if (strchr(frag, 'N') != NULL) printf("warning: N nucleotide was found in window %d of summit %lu-%lu on %s. N was not counted.\n", w+1, bs->st, bs->ed, ch1->chr); //if N is found
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
          else printf("Oops! Unexpected nucleotide %c was found in window %d of summit %lu-%lu on %s.\n", frag[i], w+1, bs->st, bs->ed, ch1->chr);
          i++;
        }

        if (!strcmp(n_flag, "A")) sprintf(tmp, "%lu\t", c_A);
        else if (!strcmp(n_flag, "T")) sprintf(tmp, "%lu\t", c_T);
        else if (!strcmp(n_flag, "C")) sprintf(tmp, "%lu\t", c_C);
        else if (!strcmp(n_flag, "G")) sprintf(tmp, "%lu\t", c_G);
        else if (!strcmp(n_flag, "AT")) sprintf(tmp, "%.3lf\t", (double)(c_A + c_T)/(double)(c_A + c_T + c_C + c_G));
        else {
          printf("error: improper n_flag: %s.\n", n_flag);
          printf("       n_flag must be either A, T, C, G or AT\n");
          goto err;
        }

        if (w == winNb - 1) tmp[strlen(tmp)-1] = '\n';
        if (strlen(ga_line_out) + strlen(tmp) + 1 > sizeof(ga_line_out)) {
          LOG("error: string per line is too long.");
          goto err;
        }
        strcat(ga_line_out, tmp); //concatenating val for each window

        if (bs->strand == '-') s = s - step; //if - strand
        else s = s + step; //if + strand or no information
      } //window

      ga_output_append (&out_head, ga_line_out);
    } //bs
  } //chromosome

  sprintf(output_name, "%s%s_nuc_hw%d_win%d_step%d_%s.txt", path, fn, hw, win, step, n_flag);

  memset(ga_line_out, '\0', LINE_STR_LEN);
  for (w = 0; w < winNb; w++) {
    sprintf(tmp, "%d\t", -hw + w*step); //header of summit file(if exist) is replaced by the relative distance header
    if (w == winNb - 1) tmp[strlen(tmp)-1] = '\n';
    if (strlen(ga_line_out) + strlen(tmp) + 1 > sizeof(ga_line_out)) {
      LOG("error: string per line is too long.");
      goto err;
    }
    strcat(ga_line_out, tmp); //concatenating relative distance
  }

  ga_write_lines (output_name, out_head, ga_line_out);

  if(ga_header_line) free(ga_header_line);
  ga_free_chr_block_fa(&chr_block_head_fa);
  ga_free_chr_block(&chr_block_head_smt);
  ga_free_chr_block(&chr_block_head_gt);
  free(frag);
  return 0;

err:
  if(ga_header_line) free(ga_header_line);
  if (chr_block_head_fa) ga_free_chr_block_fa(&chr_block_head_fa);
  if (chr_block_head_smt) ga_free_chr_block(&chr_block_head_smt);
  if (chr_block_head_gt) ga_free_chr_block(&chr_block_head_gt);
  if (frag) free(frag);
  return -1;
}

