/*
 * This program is one of the genome analysis tools.
 * This program calculates the average read distribution around summit or specific position such as TSS (Transcription Start Site).
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "sort_list.h"
#include "argument.h"
#include "ga_math.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

//static void sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, float arr[], const long smtNb, const int hw, const int step, const int win);
//static void sig_count_anti (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig_p, struct chr_block *chr_block_headsig_m, float arr[], float arr_a[], const long smtNb, const int hw, const int step, const int win);
static void sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, float arr[], const long smtNb);
static void sig_count_anti (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig_p, struct chr_block *chr_block_headsig_m, float arr[], float arr_a[], const long smtNb);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_reads_summit\n\
Summary: report the average distribution of signals around summits\n\n\
Usage:   ga_reads_summit [options] --smt <file_summit> --sig <file_signal> --sigfmt <sig format:bedgraph | sepwiggz | onewiggz>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --col_chr <int>: column number for chromosome of summit file (default:0).\n\
         --col_start <int>: column number for peak start position of summit file (default:1).\n\
         --col_end <int>: column number for peak end position of summit file (default:2).\n\
         --col_strand <int>: column number for strand of summit file (default:-1).\n\
         --header: <int> the header of summit file is preserved (default:0).\n\
         --gt: genome table file (default:NULL)\n\
         --sig_minus: signal file from minus strand. With this argument, the program calculates not only sense reads, but also anti-sense reads. (default:NULL)\n\
         --sig_d: signal denominator file like input (default:NULL)\n\
         --hw: <int> half range size (default:1000)\n\
         --step: <int> step size (default: 10)\n\
         --win: <int> window size (default:25)\n\
         --rand: <int> random simulation number. If more than 0, the simulation is performed. (default:0)\n");
  exit(0);
}

static void version()
{
  printf("'ga_reads_summit' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}


static int hf = 0;
static char *filesmt = NULL;
static char *filesig = NULL;
static char *filesig_d = NULL;
static char *filesig_m = NULL;
static char *filegenome = NULL;
static char *sigfmt = NULL;
static int col_chr = 0;
static int col_st = 1;
static int col_ed = 2;
static int col_strand = -1;
static int hw = 1000; //half window size
static int step = 10; //step size
static int win = 25; //window size
static int randnb = 0;
char *ga_header_line = NULL; //header line. Note this is global variable
static char ga_line_out[LINE_STR_LEN]; //output line including relative pos, smt_mean, CI95.00percent_U, CI95.00percent_L, smtNb, Centered, Signal

static const Argument args[] = {
  {"-h"           , ARGUMENT_TYPE_FUNCTION, usage        },
  {"--help"       , ARGUMENT_TYPE_FUNCTION, usage        },
  {"-v"           , ARGUMENT_TYPE_FUNCTION, version      },
  {"--header"     , ARGUMENT_TYPE_FLAG_ON , &hf          },
  {"--smt"        , ARGUMENT_TYPE_STRING  , &filesmt     },
  {"--sig"        , ARGUMENT_TYPE_STRING  , &filesig     },
  {"--sig_minus"  , ARGUMENT_TYPE_STRING  , &filesig_m   },
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
  {"--rand"       , ARGUMENT_TYPE_INTEGER , &randnb      },
  {NULL           , ARGUMENT_TYPE_NONE    , NULL         },
};


int main (int argc, char *argv[])
{
  argument_read(&argc, argv, args);//reading arguments
  if (filesmt == NULL || filesig == NULL || sigfmt == NULL) usage();

  struct chr_block *chr_block_headsmt = NULL; //for summit
  struct chr_block *chr_block_headsig = NULL; //for signal
  struct chr_block *chr_block_headsig_m = NULL; //for signal
  struct chr_block *chr_block_headsig_d = NULL; //for signal of denominator
  struct chr_block *chr_block_headr = NULL; //for random simulation
  struct chr_block *chr_block_headg = NULL; //for genome table

  struct chr_block *ch; //for "for loop of chr"
  struct output *output_head = NULL; //for output
  struct output *output_headr = NULL; //for output
  struct output *output_head_a = NULL; //for output
  struct output *output_headr_a = NULL; //for output

  int rel, i, r;
  float t, t2, mu_x, mu_y, ustd_y, var_x, var_y, var_xy;
  float *arr=NULL, *arr_tmp=NULL, *arr_d=NULL, *arr_tmp_d=NULL, *arr_a, *arr_tmp_a;
  float *arr_r=NULL, *arr_r_tmp=NULL, *arr_r_d=NULL, *arr_r_d_tmp=NULL, *arr_r_a=NULL, *arr_r_a_tmp=NULL;

  long smtNb, c;

  /*path, filename, and extension*/
  char path_smt[PATH_STR_LEN];
  char fn_smt[FILE_STR_LEN];
  char ext_smt[EXT_STR_LEN];
  char path_sig[PATH_STR_LEN];
  char fn_sig[FILE_STR_LEN];
  char ext_sig[EXT_STR_LEN];
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN]; //output file name

  time_t timer;


  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
Input file summit:               %s\n\
Input file signal:               %s\n\
Input file signal denominator:   %s\n\
Input file signal minus:         %s\n\
signal format:                   %s\n\
Genome file:                     %s\n\
summit col of chr, start, end:   %d, %d, %d\n\
summit col strand?:              %d\n\
half range:                      %d\n\
step size:                       %d\n\
win size:                        %d\n\
header flag:                     %d\n\
random simulation?:              %d\n\
time:                            %s\n",\
 "ga_reads_summit", filesmt, filesig, filesig_d, filesig_m, sigfmt, filegenome, col_chr, col_st, col_ed, col_strand, hw, step, win, hf, randnb, ctime(&timer) );

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
  arr_tmp = (float*)malloc(smtNb*sizeof(float)); //output arr, 1d
  if (arr_tmp == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  if (filesig_d) {//if denominator
    arr_d = (float*)malloc((((2 * hw) / step + 1) * smtNb)*sizeof(float)); //output arr, 1d
    if (arr_d == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
    arr_tmp_d = (float*)malloc(smtNb*sizeof(float)); //output arr, 1d
    if (arr_tmp_d == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
  }

  if (filesig_m) { //letting calculation of anti-strand reads mode on
    if (!strcmp(sigfmt, "bedgraph")) {
      ga_parse_bedgraph (filesig_m, &chr_block_headsig_m);
    } else if (!strcmp(sigfmt, "sepwiggz")) {
      ga_parse_sepwiggz (filesig_m, &chr_block_headsig_m);
    } else if (!strcmp(sigfmt, "onewiggz")) {
      ga_parse_onewiggz (filesig_m, &chr_block_headsig_m);
    } else {
      LOG("error: invalid signal file format.");
      goto err;
    }

    chr_block_headsig_m = ga_mergesort_chr(chr_block_headsig_m); //sorting chr of sig
    for (ch = chr_block_headsig_m; ch; ch = ch -> next) {
      ch -> sig_list = ga_mergesort_sig(ch -> sig_list); //sorting sig
    }

    //allocating arrays
    arr_a = (float*)malloc((((2 * hw) / step + 1) * smtNb)*sizeof(float)); //output arr, 1d
    if (arr_a == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
    arr_tmp_a = (float*)malloc(smtNb*sizeof(float)); //output arr, 1d
    if (arr_tmp_a == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
    sig_count_anti (chr_block_headsmt, chr_block_headsig, chr_block_headsig_m, arr, arr_a, smtNb); //counting the signal. This process is the heart of the program!
  } else {
    sig_count (chr_block_headsmt, chr_block_headsig, arr, smtNb); //counting the signal. This process is the heart of the program!
  }

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

  rel = hw; //relative pos
  t = ga_t_table (smtNb - 1); //97.5 percentile for t-dist with ddf = N -1
  t2 = t*t; //t^2

  for (i = (2 * hw) / step; i >= 0; i--) { //calculating mean, CI
    for (c = 0; c < smtNb; c++)
      arr_tmp[c] = arr[i * smtNb + c]; //choosing signal for each win by counting c = smtNb.

    mu_y = ga_mean (arr_tmp, smtNb); //mean for each win
    ustd_y = ga_ustd (arr_tmp, smtNb); //unbiased standard deviation

    if (filesig_d) { //if denominator
      for (c = 0; c < smtNb; c++)
        arr_tmp_d[c] = arr_d[i * smtNb + c]; //choosing each win
      mu_x = ga_mean (arr_tmp_d, smtNb);
      var_y = (ga_var(arr_tmp, smtNb)) / (smtNb-1);
      var_x = (ga_var(arr_tmp_d, smtNb)) / (smtNb-1);
      var_xy = (ga_covar(arr_tmp, arr_tmp_d, smtNb)) / (smtNb-1);

      if (t2 >= (mu_x * mu_x) / var_x) {
        LOG("error: normal CI cannot be calculated because denominator is not significantly different from zero.");
        goto err;
      }

      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\t%s\t%s\n", rel, mu_y / mu_x,\
         ((mu_x*mu_y - t2*var_xy)+sqrt((mu_x*mu_y - t2*var_xy)*(mu_x*mu_y - t2*var_xy)-(mu_x*mu_x - t2*var_x)*(mu_y*mu_y - t2*var_y))) / ((mu_x*mu_x) - t2*var_x),\
         ((mu_x*mu_y - t2*var_xy)-sqrt((mu_x*mu_y - t2*var_xy)*(mu_x*mu_y - t2*var_xy)-(mu_x*mu_x - t2*var_x)*(mu_y*mu_y - t2*var_y))) / ((mu_x*mu_x) - t2*var_x),\
          smtNb, fn_smt, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    }
    else { //if no denominator
      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\t%s\t%s\n", rel, mu_y, mu_y + t*ustd_y/sqrt(smtNb), mu_y - t*ustd_y/sqrt(smtNb), smtNb, fn_smt, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    }

    ga_output_add (&output_head, ga_line_out); //caution: the order is reversed
    rel -= step;
  }

  if (filesig_m) {
    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_sense.txt", path_sig, fn_sig, fn_smt, hw, win, step);
  } else {
    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d.txt", path_sig, fn_sig, fn_smt, hw, win, step);
  }
  ga_write_lines (output_name, output_head, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");

  if (filesig_m) {
    rel = hw; //relative pos

    for (i = (2 * hw) / step; i >= 0; i--) { //calculating mean, CI
      for (c = 0; c < smtNb; c++)
        arr_tmp_a[c] = arr_a[i * smtNb + c]; //choosing signal for each win by counting c = smtNb.

      mu_y = ga_mean (arr_tmp_a, smtNb); //mean for each win
      ustd_y = ga_ustd (arr_tmp_a, smtNb); //unbiased standard deviation

      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\t%s\t%s\n", rel, mu_y, mu_y + t*ustd_y/sqrt(smtNb), mu_y - t*ustd_y/sqrt(smtNb), smtNb, fn_smt, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }

      ga_output_add (&output_head_a, ga_line_out); //caution: the order is reversed
      rel -= step;
    }

    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_anti.txt", path_sig, fn_sig, fn_smt, hw, win, step);
    ga_write_lines (output_name, output_head_a, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");
  } //if (filesig_m)

  if (!randnb) { //if no random simulation, the program ends.
    goto rtfree; //goto section for free and return 0;
  }

  //the random simulation starts here.
  ga_parse_chr_bs(filegenome, &chr_block_headg, 0, 1, 1, -1, 0); //reading genome table

  arr_r = (float*)malloc((((2 * hw) / step + 1) * randnb)*sizeof(float)); //output arr, 1d
  if (arr_r == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }
  arr_r_tmp = (float*)malloc(randnb * sizeof(float)); //output arr, 1d
  if (arr_r_tmp == NULL) {
    LOG("error: lack of memory.");
    goto err;
  }

  if (filesig_d) {
    arr_r_d = (float*)malloc((((2 * hw) / step + 1) * randnb)*sizeof(float)); //output arr, 1d
    if (arr_r_d == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
    arr_r_d_tmp = (float*)malloc(randnb * sizeof(float)); //output arr, 1d
    if (arr_r_d_tmp == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
  }

  if (filesig_m) {
    arr_r_a = (float*)malloc((((2 * hw) / step + 1) * randnb)*sizeof(float)); //output arr, 1d
    if (arr_r_a == NULL) {
      LOG("error: lack of memory.");
      goto err;
    }
    arr_r_a_tmp = (float*)malloc(randnb * sizeof(float)); //output arr, 1d
    if (arr_r_a_tmp == NULL) {
      LOG("error: lack of memory.");
      goto err;
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

    if (filesig_m) {
      sig_count_anti (chr_block_headr, chr_block_headsig, chr_block_headsig_m, arr, arr_a, smtNb); //calculating signals around random postions
    } else {
      sig_count (chr_block_headr, chr_block_headsig, arr, smtNb); //calculating signals around random postions
    }
    if (filesig_d)
      sig_count (chr_block_headr, chr_block_headsig_d, arr_d, smtNb);

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

      if (filesig_m) { //if anti-sense
        for (c = 0; c < smtNb; c++)
          arr_tmp_a[c] = arr_a[i * smtNb + c]; //choosing each win
        mu_y = ga_mean (arr_tmp_a, smtNb);
        arr_r_a[i * randnb + r] = mu_y; //storing the mean
      }
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

      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\trandom\t%s\n", rel, mu_y / mu_x, mu_y / mu_x, mu_y / mu_x, smtNb, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    } else {
      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\trandom\t%s\n", rel, mu_y, mu_y, mu_y, smtNb, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
    }
    ga_output_add (&output_headr, ga_line_out); //caution: the order is reversed

    if (filesig_m) {
      for (r = 0; r < randnb; r++)
        arr_r_a_tmp[r] = arr_r_a[i * randnb + r]; //choosing signal for each win by counting r = randnb.
      mu_y = ga_mean (arr_r_a_tmp, randnb); //mean for each win

      if (sprintf(ga_line_out, "%d\t%f\t%f\t%f\t%ld\trandom\t%s\n", rel, mu_y, mu_y, mu_y, smtNb, fn_sig) == EOF) {
        LOG("error: the summit name or signal name is too long.");
        goto err;
      }
      ga_output_add (&output_headr_a, ga_line_out); //caution: the order is reversed
    }

    rel -= step;
  } //i

  if (filesig_m) {
    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_sense_random%d.txt", path_sig, fn_sig, fn_smt, hw, win, step, randnb);
    ga_write_lines (output_name, output_headr, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");
    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_anti_random%d.txt", path_sig, fn_sig, fn_smt, hw, win, step, randnb);
    ga_write_lines (output_name, output_headr_a, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");
  } else {
    sprintf(output_name, "%s%s_around_%s_halfwid%dwinsize%dstep%d_random%d.txt", path_sig, fn_sig, fn_smt, hw, win, step, randnb);
    ga_write_lines (output_name, output_headr, "relative_pos\tsmt_mean\tCI95.00percent_U\tCI95.00percent_L\tsmtNb\tCentered\tSignal\n");
  }

  goto rtfree;

rtfree:
  if (arr) free(arr);
  if (arr_tmp) free(arr_tmp);
  if (arr_d) free(arr_d);
  if (arr_tmp_d) free(arr_tmp_d);
  if (arr_a) free(arr_a);
  if (arr_tmp_a) free(arr_tmp_a);
  if (arr_r) free(arr_r);
  if (arr_r_tmp) free(arr_r_tmp);
  if (arr_r_d) free(arr_r_d);
  if (arr_r_d_tmp) free(arr_r_d_tmp);
  if (arr_r_a) free(arr_r_a);
  if (arr_r_a_tmp) free(arr_r_a_tmp);
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (chr_block_headsig_m) ga_free_chr_block(&chr_block_headsig_m);
  if (chr_block_headr) ga_free_chr_block(&chr_block_headr);
  if (chr_block_headg) ga_free_chr_block(&chr_block_headg);
  if (output_head) ga_free_output(&output_head);
  if (output_headr) ga_free_output(&output_headr);
  if (output_head_a) ga_free_output(&output_head_a);
  if (output_headr_a) ga_free_output(&output_headr_a);
  if (ga_header_line) free(ga_header_line);

  return 0;

err:
  if (arr) free(arr);
  if (arr_tmp) free(arr_tmp);
  if (arr_d) free(arr_d);
  if (arr_tmp_d) free(arr_tmp_d);
  if (arr_a) free(arr_a);
  if (arr_tmp_a) free(arr_tmp_a);
  if (arr_r) free(arr_r);
  if (arr_r_tmp) free(arr_r_tmp);
  if (arr_r_d) free(arr_r_d);
  if (arr_r_d_tmp) free(arr_r_d_tmp);
  if (arr_r_a) free(arr_r_a);
  if (arr_r_a_tmp) free(arr_r_a_tmp);
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (chr_block_headsig_m) ga_free_chr_block(&chr_block_headsig_m);
  if (chr_block_headr) ga_free_chr_block(&chr_block_headr);
  if (chr_block_headg) ga_free_chr_block(&chr_block_headg);
  if (output_head) ga_free_output(&output_head);
  if (output_headr) ga_free_output(&output_headr);
  if (output_head_a) ga_free_output(&output_head_a);
  if (output_headr_a) ga_free_output(&output_headr_a);
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

        if (j1_tmp == NULL && bs->strand != '-') { //if j1_tmp is not set for the chr
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


static void sig_count_anti (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig_p, struct chr_block *chr_block_headsig_m, float arr[], float arr_a[], const long smtNb)
{
  struct chr_block *ch_smt, *ch_sig;
  struct bs *bs;
  struct sig *j1, *j1_tmp = NULL; //j1 is the pointer to chr_block_headsig which is counted in the window. j1_tmp is the 'memory' of j1 which act as the marker of the previous position of j1 to speed up the calculation. Thanks to j1_tmp, we don't have to search the signal position of 1 for each chr, rather we can start the searching from the previous position.
  int i, fl, winNb = (2 * hw) / step + 1;
  long c=0, c_tmp=0, st, ed, tmp_st, tmp_ed;
  float val_tmp;

  for (ch_smt = chr_block_headsmt; ch_smt; ch_smt = ch_smt->next) {
//    printf("calculating reads on %s\n", ch_smt->chr);
    for (ch_sig = chr_block_headsig_p; ch_sig; ch_sig = ch_sig->next) {
      if (!strcmp(ch_smt->chr, ch_sig->chr)) break; //if the same chr is included in smt and sig
    }

    if (ch_sig == NULL) { //if chr in smt is not included in sig...
      for (bs = ch_smt->bs_list; bs; bs = bs->next) {
        for (i = 0; i < winNb; i++) {
          arr[i * smtNb + c] = 0.0; //assigning value 0.0 if chr in smt is not included in sig. 
          arr_a[i * smtNb + c] = 0.0; //assigning value 0.0 if chr in smt is not included in sig. 
        }
        c++; //counting up for each bs
      }
      continue;
    }

    c_tmp = c; //memory of c

    //calculating for plus strand reads
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
          if (bs->strand == '-') arr_a[(winNb -1 - i) * smtNb + c] = 0.0; //assigning value 0.0
          else arr[i * smtNb + c] = 0.0; //assigning value 0.0
          st += step;
          ed += step;
          continue;
        }

        if (j1_tmp == NULL && bs->strand != '-') { //if j1_tmp is not set for the chr
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
        if (bs->strand == '-') arr_a[(winNb -1 - i) * smtNb + c] = val_tmp / (float)win;
        else arr[i * smtNb + c] = val_tmp / (float)win;
        st += step;
        ed += step;
      }
      c++; //counting up for each bs
    } //bs

    c = c_tmp; //memory of c
    for (ch_sig = chr_block_headsig_m; ch_sig; ch_sig = ch_sig->next) {
      if (!strcmp(ch_smt->chr, ch_sig->chr)) break; //if the same chr is included in smt and sig
    }

    //calculating for minus strand reads
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
          else arr_a[i * smtNb + c] = 0.0; //assigning value 0.0
          st += step;
          ed += step;
          continue;
        }

        if (j1_tmp == NULL && bs->strand != '-') { //if j1_tmp is not set for the chr
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
        else arr_a[i * smtNb + c] = val_tmp / (float)win;
        st += step;
        ed += step;
      }
      c++; //counting up for each bs
    } //bs

  } //chr

  return;
}



