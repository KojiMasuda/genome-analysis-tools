/*
 * This program is one of the genome analysis tools.
 * This program calculates read distributions around ALL summits or specific positions such as TSS (Transcription Start Site).
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

//static void sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, float arr[], const long smtNb, const int hw, const int step, const int win);
static void sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, float arr[], const long smtNb);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_reads_summit_all\n\
Summary: report the distribution of all signals around summits\n\n\
Usage:   ga_reads_summit_all [options] --smt <file_summit> --sig <file_signal> --sigfmt <sig format:bedgraph | sepwiggz | onewiggz>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --col_chr <int>: column number for chromosome of summit file (default:0).\n\
         --col_start <int>: column number for peak start position of summit file (default:1).\n\
         --col_end <int>: column number for peak end position of summit file (default:2).\n\
         --col_strand <int>: column number for strand of summit file (default:-1).\n\
         --header: the first line of summit file is considered as header (default:off).\n\
         --gt: genome table file (default:NULL)\n\
         --sig_d: signal denominator file like input (default:NULL)\n\
         --hw: <int> half range size (default:1000)\n\
         --step: <int> step size (default: 10)\n\
         --win: <int> window size (default:25)\n");
  exit(0);
}

static void version()
{
  printf("'ga_reads_summit_all' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}


static int hf = 0;
static char *filesmt = NULL;
static char *filesig = NULL;
static char *filesig_d = NULL;
static char *filegenome = NULL;
static char *sigfmt = NULL;
static int col_chr = 0;
static int col_st = 1;
static int col_ed = 2;
static int col_strand = -1;
static int hw = 1000; //half window size
static int step = 10; //step size
static int win = 25; //window size
char *ga_header_line = NULL; //header line. Note this is external global variable
static char ga_line_out[LINE_STR_LEN]; //output line including relative pos, smt_mean, CI95.00percent_U, CI95.00percent_L, smtNb, Centered, Signal

static const Argument args[] = {
  {"-h"           , ARGUMENT_TYPE_FUNCTION, usage        },
  {"--help"       , ARGUMENT_TYPE_FUNCTION, usage        },
  {"-v"           , ARGUMENT_TYPE_FUNCTION, version      },
  {"--header"     , ARGUMENT_TYPE_FLAG_ON , &hf          },
  {"--smt"        , ARGUMENT_TYPE_STRING  , &filesmt     },
  {"--sig"        , ARGUMENT_TYPE_STRING  , &filesig     },
  {"--sig_d"      , ARGUMENT_TYPE_STRING  , &filesig_d   },
  {"--gt"         , ARGUMENT_TYPE_STRING  , &filegenome  },
  {"--sigfmt"     , ARGUMENT_TYPE_STRING  , &sigfmt      },
  {"--col_chr"    , ARGUMENT_TYPE_INTEGER , &col_chr     },
  {"--col_start"  , ARGUMENT_TYPE_INTEGER , &col_st      },
  {"--col_end"    , ARGUMENT_TYPE_INTEGER , &col_ed      },
  {"--col_strand" , ARGUMENT_TYPE_INTEGER , &col_strand  },
  {"--hw"         , ARGUMENT_TYPE_INTEGER , &hw          },
  {"--step"       , ARGUMENT_TYPE_INTEGER , &step        },
  {"--win"        , ARGUMENT_TYPE_INTEGER , &win         },
  {NULL           , ARGUMENT_TYPE_NONE    , NULL         },
};


int main (int argc, char *argv[])
{
  argument_read(&argc, argv, args);//reading arguments
  if (filesmt == NULL || filesig == NULL || sigfmt == NULL) usage();

  struct chr_block *chr_block_headsmt = NULL; //for summit
  struct chr_block *chr_block_headsig = NULL; //for signal
  struct chr_block *chr_block_headsig_d = NULL; //for signal of denominator

  struct chr_block *ch; //for "for loop of chr"
  struct output *output_head = NULL; //for output

  int rel, i;
  float *arr=NULL, *arr_d=NULL; //, *arr_a, *arr_tmp_a;

  long smtNb, c;

  /*path, filename, and extension*/
  char path_smt[PATH_STR_LEN];
  char fn_smt[FILE_STR_LEN];
  char ext_smt[EXT_STR_LEN];
  char path_sig[PATH_STR_LEN];
  char fn_sig[FILE_STR_LEN];
  char ext_sig[EXT_STR_LEN];
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN]; //output file name
  char str_tmp[32]; //for each value with \t

  time_t timer;


  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
Input file summit:               %s\n\
Input file signal:               %s\n\
Input file signal denominator:   %s\n\
signal format:                   %s\n\
summit col of chr, start, end:   %d, %d, %d\n\
summit col strand?:              %d\n\
half range:                      %d\n\
step size:                       %d\n\
win size:                        %d\n\
header flag:                     %d\n\
time:                            %s\n",\
 "ga_reads_summit_all", filesmt, filesig, filesig_d, sigfmt, col_chr, col_st, col_ed, col_strand, hw, step, win, hf, ctime(&timer) );

  ga_parse_file_path (filesmt, path_smt, fn_smt, ext_smt); //parsing input file name into path, file name, and extension
  ga_parse_file_path (filesig, path_sig, fn_sig, ext_sig);

  ga_parse_chr_bs(filesmt, &chr_block_headsmt, col_chr, col_st, col_ed, col_strand, hf); //parsing each binding sites for each chromosome

  // reading signal file
  if (!strcmp(sigfmt, "bedgraph")) {
    ga_parse_bedgraph (filesig, &chr_block_headsig);
  } else if (!strcmp(sigfmt, "sepwiggz")) {
    ga_parse_sepwiggz (filesig, &chr_block_headsig);
  } else if (!strcmp(sigfmt, "onewiggz")) {
    ga_parse_onewiggz (filesig, &chr_block_headsig);
  } else {
    LOG("error: invalid signal file format.");
    goto err;
  }

  // sorting summit and sig
  chr_block_headsmt = ga_mergesort_chr(chr_block_headsmt);
  chr_block_headsig = ga_mergesort_chr(chr_block_headsig);

  for (ch = chr_block_headsmt; ch; ch = ch -> next) {
    ch -> bs_list = ga_mergesort_bs(ch -> bs_list);
  }
  for (ch = chr_block_headsig; ch; ch = ch -> next) {
    ch -> sig_list = ga_mergesort_sig(ch -> sig_list);
  }


  smtNb = ga_count_peaks (chr_block_headsmt); //counting smt number
  printf("smtnb:%ld\n", smtNb);

  //allocating arrays
  arr = (float*)malloc((((2 * hw) / step + 1) * smtNb)*sizeof(float)); //output arr, 1d
  if (arr == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  if (filesig_d) {//if denominator
    arr_d = (float*)malloc((((2 * hw) / step + 1) * smtNb)*sizeof(float)); //output arr, 1d
    if (arr_d == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
  }

  sig_count (chr_block_headsmt, chr_block_headsig, arr, smtNb); //counting the signal. This process is the heart of the program!

  if (filesig_d) {//if denominator
    if (!strcmp(sigfmt, "bedgraph")) {
      ga_parse_bedgraph (filesig_d, &chr_block_headsig_d);
    } else if (!strcmp(sigfmt, "sepwiggz")) {
      ga_parse_sepwiggz (filesig_d, &chr_block_headsig_d);
    } else if (!strcmp(sigfmt, "onewiggz")) {
      ga_parse_onewiggz (filesig_d, &chr_block_headsig_d);
    } else {
      LOG("error: invalid signal file format.");
      goto err;
    }

    chr_block_headsig_d = ga_mergesort_chr(chr_block_headsig_d); //sorting chr of sig
    for (ch = chr_block_headsig_d; ch; ch = ch -> next) {
      ch -> sig_list = ga_mergesort_sig(ch -> sig_list); //sorting sig
    }
    sig_count (chr_block_headsmt, chr_block_headsig_d, arr_d, smtNb);
  }

  if (filesig_d) {
    for (c = smtNb - 1; c >= 0 ; c--) { 
      memset(ga_line_out, '\0', sizeof(ga_line_out)); //assigning \0 into ga_line_out

      for (i = 0 ;i <= (2 * hw) / step ;i++) { //concatenating val for each window
        if (arr_d[i * smtNb + c] == 0 && i == ((2 * hw) / step)) { //if denominator is zero...
          printf("warning: denominator of win%d for peak%ld was zero.\n", i+1, c+1);
          sprintf(str_tmp, "NA\n");
        } else if (arr_d[i * smtNb + c] == 0) { //if denominator is zero...
          printf("warning: denominator of win%d for peak%ld was zero.\n", i+1, c+1);
          sprintf(str_tmp, "NA\t");
        } else if (i == ((2 * hw) / step)) {
          sprintf(str_tmp, "%f\n", arr[i * smtNb + c] / arr_d[i * smtNb + c]); //adding val for the last window
        } else {
          sprintf(str_tmp, "%f\t", arr[i * smtNb + c] / arr_d[i * smtNb + c]); //adding val
        }
        if (strlen(ga_line_out) + strlen(str_tmp) + 1 > sizeof(ga_line_out)) {
          LOG("error: string per line is too long.");
          goto err;
        }
        strcat(ga_line_out, str_tmp); //concatenating val for each window
      }
      ga_output_add (&output_head, ga_line_out); //caution: the order is reversed
    }
  } else { //if you want the peak ID, you have to locate the peak by accessing chr_block_headsmt, invert c, and use ga_output_append...
    for (c = smtNb - 1; c >= 0 ; c--) {
      memset(ga_line_out, '\0', sizeof(ga_line_out)); //assigning \0 into ga_line_out

      for (i = 0 ;i <= (2 * hw) / step ;i++) { //concatenating val for each window
        if (i == ((2 * hw) / step)) {
          sprintf(str_tmp, "%f\n", arr[i * smtNb + c]); //adding val for the last window
        } else {
          sprintf(str_tmp, "%f\t", arr[i * smtNb + c]); //adding val
        }
        if (strlen(ga_line_out) + strlen(str_tmp) + 1 > sizeof(ga_line_out)) {
          LOG("error: string per line is too long.");
          goto err;
        }
        strcat(ga_line_out, str_tmp); //concatenating val for each window
      }
      ga_output_add (&output_head, ga_line_out); //caution: the order is reversed
    }
  }

  sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_all.txt", path_sig, fn_sig, fn_smt, hw, win, step);
  rel = -hw; //relative pos
  memset(ga_line_out, '\0', sizeof(ga_line_out)); //assigning \0 into ga_line_out

  for (i = 0 ;i <= (2 * hw) / step ;i++) { //concatenating relative positions for each window
    if (i == ((2 * hw) / step)) {
      sprintf(str_tmp, "%d\n", rel); //relative position for the last window
    } else {
      sprintf(str_tmp, "%d\t", rel); //relative position
    }
    if (strlen(ga_line_out) + strlen(str_tmp) + 1 > sizeof(ga_line_out)) {
      LOG("error: string per line is too long.");
      goto err;
    }
    strcat(ga_line_out, str_tmp); //concatenating val for each window
    rel += step;
  }
  ga_write_lines (output_name, output_head, ga_line_out);

/* TEMP!
  if (!randnb) { //if no random simulation, the program ends.
    goto rtfree; //goto section for free and return 0;
  }

  //the random simulation starts here.
  ga_parse_chr_bs(filegenome, &chr_block_headg, 0, 1, 1, -1, 0); //reading genome table

  arr_r = (float*)malloc((((2 * hw) / step + 1) * randnb)*sizeof(float)); //output arr, 1d
  if (arr_r == NULL) {
    LOG("error: lack of memory.");
  }
  arr_r_tmp = (float*)malloc(randnb * sizeof(float)); //output arr, 1d
  if (arr_r_tmp == NULL) {
    LOG("error: lack of memory.");
  }

  if (filesig_d) {
    arr_r_d = (float*)malloc((((2 * hw) / step + 1) * randnb)*sizeof(float)); //output arr, 1d
    if (arr_r_d == NULL) {
      LOG("error: lack of memory.");
    }
    arr_r_d_tmp = (float*)malloc(randnb * sizeof(float)); //output arr, 1d
    if (arr_r_d_tmp == NULL) {
      LOG("error: lack of memory.");
    }
  }

  for (r = 0; r < randnb; r++) {
    printf("\rsimulation cycle: %d", r + 1);
    fflush(stdout);
    if (r) { //if chr_block_headr is alreadly allocated, it should be freed.
      ga_free_chr_block(&chr_block_headr); //free chr_block_smt
      chr_block_headr = NULL; //
    }

    ga_parse_chr_bs_rand (&chr_block_headr, chr_block_headsmt, chr_block_headg, hw); //picking up random positions

    chr_block_headr = ga_mergesort_chr(chr_block_headr);
    for (ch = chr_block_headr; ch; ch = ch -> next) {
      ch -> bs_list = ga_mergesort_bs(ch -> bs_list); //sorting
    }

    sig_count (chr_block_headr, chr_block_headsig, arr, smtNb, hw, step, win); //calculating signals around random postions
//    }
    if (filesig_d)
      sig_count (chr_block_headr, chr_block_headsig_d, arr_d, smtNb, hw, step, win);

    for (i = 0; i < (2 * hw) / step + 1; i++) { //calculating mean
      for (c = 0; c < smtNb; c++)
        arr_tmp[c] = arr[i * smtNb + c]; //choosing signal for each win by counting c = smtNb.
      mu_y = ga_mean (arr_tmp, smtNb); //mean for each win

      if (filesig_d) { //if denominator
        for (c = 0; c < smtNb; c++)
          arr_tmp_d[c] = arr_d[i * smtNb + c]; //choosing each win
        mu_x = ga_mean (arr_tmp_d, smtNb);
      }


      arr_r[i * randnb + r] = mu_y; //storing the mean
      if (filesig_d)
        arr_r_d[i * randnb + r] = mu_x; //storing the mean

    } //i
  } //r
  printf("\n");

  rel = hw; //relative pos
  for (i = (2 * hw) / step; i >= 0; i--) {
    for (r = 0; r < randnb; r++)
      arr_r_tmp[r] = arr_r[i * randnb + r]; //choosing signal for each win by counting r = randnb.
    mu_y = ga_mean (arr_r_tmp, randnb); //mean for each win

    if (filesig_d) { //if denominator
      for (r = 0; r < randnb; r++)
        arr_r_d_tmp[r] = arr_r_d[i * randnb + r]; //choosing signal for each win by counting r = randnb.
      mu_x = ga_mean (arr_r_d_tmp, randnb); //mean for each win

      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%lu\trandom\t%s\n", rel, mu_y / mu_x, mu_y / mu_x, mu_y / mu_x, smtNb, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    } else {
      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%lu\trandom\t%s\n", rel, mu_y, mu_y, mu_y, smtNb, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    }
    ga_output_add (&output_headr, ga_line_out);


    rel -= step;
  } //i

  sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_random%d.txt", path_sig, fn_sig, fn_smt, hw, win, step, randnb);
  ga_write_lines (output_name, output_headr, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");
//  }
*/

  goto rtfree;

rtfree:
  if (arr) free(arr);
  if (arr_d) free(arr_d);
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);

  return 0;

err:
  if (arr) free(arr);
  if (arr_d) free(arr_d);
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);
  return -1;
}

//the structure of arr is [win1:peak1,peak2...peakN|win2:peak1,peak2...peakN|...|winN:peak1,peak2...peakN]
//the structure of arr_r is [r1:peak1,peak2...peakN|r2:peak1,peak2...peakN|...|rN:peak1,peak2...peakN]
static void sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, float arr[], const long smtNb)
{
  struct chr_block *ch_smt, *ch_sig;
  struct bs *bs;
  struct sig *j1, *j1_tmp = NULL; //j1 is the pointer to chr_block_headsig which is counted in the window. j1_tmp is the 'memory' of j1 which act as the marker of the previous position of j1 to speed up the calculation. Thanks to j1_tmp, we don't have to search the signal position of 1 for each chr, rather we can start the searching from the previous position.
  int i, fl, winNb = (2 * hw) / step + 1;
  long c=0, st, ed, tmp_st, tmp_ed;
  float val_tmp;

  for (ch_smt = chr_block_headsmt; ch_smt; ch_smt = ch_smt->next) {
//    printf("calculating reads on %s\n", ch_smt->chr);
    for (ch_sig = chr_block_headsig; ch_sig; ch_sig = ch_sig->next) {
      if (!strcmp(ch_smt->chr, ch_sig->chr)) break; //if the same chr is included in smt and sig
    }

    if (ch_sig == NULL) { //if chr in smt is not included in sig...
      for (bs = ch_smt->bs_list; bs; bs = bs->next) {
        for (i = 0; i < winNb; i++) arr[i * smtNb + c] = 0.0; //assigning value 0.0 if chr in smt is not included in sig. 
        c++; //counting up for each bs
      }
      continue;
    }

    j1_tmp = NULL; //the "marker" of signal position to speed up the calc. j1_tmp is the left most position for each bs.
    for (bs = ch_smt->bs_list; bs; bs = bs->next) {
      fl = 0; //init of flag
      if (bs->strand == '-') {//if the summit is on minus strand
        st = bs->ed - hw - win / 2; //start pos
        ed = bs->ed - hw + win / 2; //end pos
      } else {
        st = bs->st - hw - win / 2; //start pos
        ed = bs->st - hw + win / 2; //end pos
      }
      for (i = 0; i < winNb; i++) {
        if (j1_tmp == NULL) {
          for (j1 = ch_sig->sig_list; j1; j1=j1->next) { //here's the slowest part...
            if (st < j1->ed && j1->st < ed) {
              break; //if one of sig block is inside the win
            } else if (j1->st >= ed) { //if there's no chance for j1 to overlap win
              j1 = NULL;
              break;
            }
          }
        }
        else {
          for (j1 = j1_tmp; j1; j1=j1->next) { //here's the slowest part...
            if (st < j1->ed && j1->st < ed) {
              break; //if one of sig block is inside the win
            } else if (j1->st >= ed) { //if there's no chance for j1 to overlap win
              j1 = NULL;
              break;
            }
          }
        }

        if (j1 == NULL) { //if the win is the right side of the most right sig block
          if (bs->strand == '-') arr[(winNb -1 - i) * smtNb + c] = 0.0; //assigning value 0.0
          else arr[i * smtNb + c] = 0.0; //assigning value 0.0
          st += step;
          ed += step;
          continue;
        }

        if (j1_tmp == NULL && bs->strand != '-') { //if j1_tmp is not set for the chr and the strand is not minus.
          j1_tmp = j1;
          fl = 1;
        }
        else if (!fl && bs->strand != '-') { //if j1_tmp is not set for the bs (fl == 0) and the strand is not minus.
          j1_tmp = j1;
          fl = 1;
        }
        
        val_tmp = 0;
        for (;j1 ; j1 = j1->next) {
          if (j1->st >= ed) break; //if the sig pos is out of the win
          if (st > j1->st) tmp_st = st; //if st of sig block is up-stream pos of st
          else tmp_st = j1->st;
          if (j1->ed > ed) tmp_ed = ed; //if ed of sig block is down-stream pos of ed
          else tmp_ed = j1->ed;
          val_tmp += (j1->val) * (tmp_ed - tmp_st); //adding the val*len of sig block
        }
        if (bs->strand == '-') arr[(winNb -1 - i) * smtNb + c] = val_tmp / (float)win;
        else arr[i * smtNb + c] = val_tmp / (float)win;
        st += step;
        ed += step;
      }
      c++; //counting up for each bs
    }
  }

  return;
}

