/*
 * This program is one of the genome analysis tools.
 * This program calculates signals around summits(with fixed length) or regions.
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

//static int sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, struct chr_block *chr_block_headsig_d, struct output **output_head, const int hw, const char *region_mode);
static int sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, struct chr_block *chr_block_headsig_d, struct output **output_head);

static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_reads_region\n\
Summary: report the signals around summits or regions\n\n\
Usage:   ga_reads_region [options] --smt <file_summit> --sig <file_signal> --sigfmt <sig format:bedgraph | sepwiggz | onewiggz> --mode smt --col_smt <column of summit>\n\
   or:   ga_reads_region [options] --smt <file_summit> --sig <file_signal> --sigfmt <sig format:bedgraph | sepwiggz | onewiggz> --mode <region mode: region | up-tss | tss-dw | up-tss-dw | up-tes | tes-dw | up-tes-dw>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --col_chr <int>: column number for chromosome of summit file (default:0).\n\
         --col_start <int>: column number for peak start position of summit file (default:1).\n\
         --col_end <int>: column number for peak end position of summit file (default:2).\n\
         --col_smt <int>: column number for summit position of summit file (default:1).\n\
         --col_strand <int>: column number for strand of summit file (default:-1).\n\
         --header: <int> the header of summit file is preserved (default:0).\n\
         --sig_d: signal denominator file like input (default:NULL)\n\
         --hw: <int> half range size (default:1000)\n");
  exit(0);
}

static void version()
{
  printf("'ga_reads_region' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}


static int hf = 0;
static char *filesmt = NULL;
static char *filesig = NULL;
static char *filesig_d = NULL;
static char *sigfmt = NULL;
static char *region_mode = NULL;
static int col_chr = 0;
static int col_st = 1;
static int col_ed = 2;
static int col_strand = -1;
static int hw = 1000; //half window size
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
  {"--sigfmt"     , ARGUMENT_TYPE_STRING  , &sigfmt      },
  {"--mode"       , ARGUMENT_TYPE_STRING  , &region_mode },
  {"--col_chr"    , ARGUMENT_TYPE_INTEGER , &col_chr     },
  {"--col_start"  , ARGUMENT_TYPE_INTEGER , &col_st      },
  {"--col_end"    , ARGUMENT_TYPE_INTEGER , &col_ed      },
  {"--col_strand" , ARGUMENT_TYPE_INTEGER , &col_strand  },
  {"--hw"         , ARGUMENT_TYPE_INTEGER , &hw          },
  {"--col_smt"    , ARGUMENT_TYPE_INTEGER , &col_st      },
  {NULL           , ARGUMENT_TYPE_NONE    , NULL         },
};


int main (int argc, char *argv[])
{
  argument_read(&argc, argv, args);//reading arguments
  if (filesmt == NULL || filesig == NULL || sigfmt == NULL || region_mode == NULL) usage();

  struct chr_block *chr_block_headsmt = NULL; //for summit
  struct chr_block *chr_block_headsig = NULL; //for signal
  struct chr_block *chr_block_headsig_d = NULL; //for signal of denominator

  struct chr_block *ch; //for "for loop of chr"
  struct output *output_head = NULL; //for output

  /*path, filename, and extension*/
  char path_smt[PATH_STR_LEN];
  char fn_smt[FILE_STR_LEN];
  char ext_smt[EXT_STR_LEN];
  char path_sig[PATH_STR_LEN];
  char fn_sig[FILE_STR_LEN];
  char ext_sig[EXT_STR_LEN];
  char path_sig_d[PATH_STR_LEN];
  char fn_sig_d[FILE_STR_LEN];
  char ext_sig_d[EXT_STR_LEN];
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN]; //output file name

  time_t timer;


  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
Input file summit:               %s\n\
Input file signal:               %s\n\
Input file signal denominator:   %s\n\
signal format:                   %s\n\
region mode:                     %s\n\
summit col of chr, start, end:   %d, %d, %d\n\
summit col of summit:            %d\n\
summit col strand?:              %d\n\
half range:                      %d\n\
header flag:                     %d\n\
time:                            %s\n",\
 "ga_reads_region", filesmt, filesig, filesig_d, sigfmt, region_mode, col_chr, col_st, col_ed, col_st, col_strand, hw, hf, ctime(&timer) );

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

  if (filesig_d) {//if denominator
    ga_parse_file_path (filesig_d, path_sig_d, fn_sig_d, ext_sig_d);

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
  }

  if (sig_count (chr_block_headsmt, chr_block_headsig, chr_block_headsig_d, &output_head) != 0) {
    LOG("error: in sig_count function.");
    goto err;
  }

  if (filesig_d) { //if denominator
    sprintf(output_name, "%s%s_around_%s_halfwid%d_mode_%s_devided%s.txt", path_sig, fn_sig, fn_smt, hw, region_mode, fn_sig_d);
  } else {
    sprintf(output_name, "%s%s_around_%s_halfwid%d_mode_%s.txt", path_sig, fn_sig, fn_smt, hw, region_mode);
  }

  if (ga_header_line != NULL) { //if header line
    if (add_one_val(ga_line_out, ga_header_line, "signal_region\n") < 0) {
      LOG("error: output line was too long.");
      goto err; //adding one extra column
    }
    ga_write_lines (output_name, output_head, ga_line_out); //note that header is line_out, not ga_header_line
  }
  else ga_write_lines (output_name, output_head, ga_header_line);

  goto rtfree;

rtfree:
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);

  return 0;

err:
  if (chr_block_headsmt) ga_free_chr_block(&chr_block_headsmt);
  if (chr_block_headsig) ga_free_chr_block(&chr_block_headsig);
  if (chr_block_headsig_d) ga_free_chr_block(&chr_block_headsig_d);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) free(ga_header_line);
  return -1;
}

//static int sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, struct chr_block *chr_block_headsig_d, struct output **output_head, const int hw, const char *region_mode)
static int sig_count (struct chr_block *chr_block_headsmt, struct chr_block *chr_block_headsig, struct chr_block *chr_block_headsig_d, struct output **output_head)
{
  struct chr_block *ch_smt, *ch_sig, *ch_sig_d = NULL;
  struct bs *bs;
  struct sig *j1, *j1_tmp = NULL, *j1_d, *j1_tmp_d = NULL; //j1 is the pointer to chr_block_headsig which is counted in the window. j1_tmp is the 'memory' of j1 which act as the marker of the previous position of j1 to speed up the calculation. Thanks to j1_tmp, we don't have to search the signal position of 1 for each chr, rather we can start the searching from the previous position.
  int fl, fl_d;
  long st, ed, tmp_st, tmp_ed;
  float val_tmp, val_tmp_d;
  char tmp[64], tag[10];

  for (ch_smt = chr_block_headsmt; ch_smt; ch_smt = ch_smt->next) {
    for (ch_sig = chr_block_headsig; ch_sig; ch_sig = ch_sig->next) {
      if (!strcmp(ch_smt->chr, ch_sig->chr)) break; //if the same chr is included in smt and sig
    }
    if (chr_block_headsig_d) {
      for (ch_sig_d = chr_block_headsig_d; ch_sig_d; ch_sig_d = ch_sig_d->next) {
        if (!strcmp(ch_smt->chr, ch_sig_d->chr)) break; //if the same chr is included in smt and sig
      }
    }

    if (ch_sig == NULL) { //if chr in smt is not included in sig...
      for (bs = ch_smt->bs_list; bs; bs = bs->next) {
        sprintf(tmp, "%f\n", 0.0);
        if (add_one_val(ga_line_out, bs->line, tmp) != 0){
          LOG("error: output line was too long.");
          return -1;
        }
        ga_output_append (output_head, ga_line_out);
      }
      continue;
    } else if (chr_block_headsig_d && ch_sig_d == NULL) {
      printf("error: %s is not found in signal denominator.\n", ch_smt->chr);
      return -1;
    }

    j1_tmp = NULL; //the "marker" of signal position to speed up the calc. j1_tmp is the left most position for each bs.
    j1_tmp_d = NULL; //the "marker" of signal position to speed up the calc. j1_tmp is the left most position for each bs.
    for (bs = ch_smt->bs_list; bs; bs = bs->next) {
      fl = 0; //init of flag
      fl_d = 0; //init of flag

      if (!strcmp(region_mode, "smt")) { //if summit +- fixed half window
        st = bs->st - hw; //start pos
        ed = bs->st + hw; //end pos
        sprintf(tag, "tss"); //temporary tag...
      } else if (!strcmp(region_mode, "region")) {
        st = bs->st; //start pos
        ed = bs->ed; //end pos
        sprintf(tag, "tss"); //temporary tag...
      } else if (!strcmp(region_mode, "up-tss")) {
        if (bs->strand == '-') {//if the region is on minus strand
          st = bs->ed; //start pos
          ed = bs->ed + hw * 2; //end pos
        } else {
          st = bs->st - hw * 2; //start pos
          ed = bs->st; //end pos
        }
        sprintf(tag, "tss"); //tag
      } else if (!strcmp(region_mode, "tss-dw")) {
        if (bs->strand == '-') {//if the region is on minus strand
          st = bs->ed - hw * 2; //start pos
          ed = bs->ed; //end pos
        } else {
          st = bs->st; //start pos
          ed = bs->st + hw * 2; //end pos
        }
        sprintf(tag, "tss"); //tag
      } else if (!strcmp(region_mode, "up-tss-dw")) {
        if (bs->strand == '-') {//if the region is on minus strand
          st = bs->ed - hw * 2; //start pos
          ed = bs->ed + hw * 2; //end pos
        } else {
          st = bs->st - hw * 2; //start pos
          ed = bs->st + hw * 2; //end pos
        }
        sprintf(tag, "tss"); //tag
      } else if (!strcmp(region_mode, "up-tes")) {
        if (bs->strand == '+') {//if the region is on minus strand
          st = bs->ed - hw * 2; //start pos
          ed = bs->ed; //end pos
        } else {
          st = bs->st; //start pos
          ed = bs->st + hw * 2; //end pos
        }
        sprintf(tag, "tes"); //tag
      } else if (!strcmp(region_mode, "tes-dw")) {
        if (bs->strand == '+') {//if the region is on minus strand
          st = bs->ed; //start pos
          ed = bs->ed + hw * 2; //end pos
        } else {
          st = bs->st - hw * 2; //start pos
          ed = bs->st; //end pos
        }
        sprintf(tag, "tes"); //tag
      } else if (!strcmp(region_mode, "up-tes-dw")) {
        if (bs->strand == '+') {//if the region is on minus strand
          st = bs->ed - hw * 2; //start pos
          ed = bs->ed + hw * 2; //end pos
        } else {
          st = bs->st - hw * 2; //start pos
          ed = bs->st + hw * 2; //end pos
        }
        sprintf(tag, "tes"); //tag
      } else {
        printf("error: input correct region_mode. Your mode: %s\n", region_mode);
        return -1;
      }

      if (j1_tmp == NULL) j1 = ch_sig->sig_list;
      else j1 = j1_tmp;

      for (; j1; j1=j1->next) { //here's the slowest part...
        if (st < j1->ed && j1->st < ed) {
          break; //if one of sig block is inside the win
        } else if (j1->st >= ed) { //if there's no chance for j1 to overlap win
          j1 = NULL;
          break;
        }
      }

      if (j1 == NULL) { //if the win is the right side of the most right sig block
        sprintf(tmp, "%f\n", 0.0);
        if (add_one_val(ga_line_out, bs->line, tmp) != 0){
          LOG("error: output line was too long.");
          return -1; //making output link list with NA
        }
        ga_output_append (output_head, ga_line_out);
        continue;
      }

      if (!fl && !strcmp(tag, "tss") && bs->strand != '-') { //if j1_tmp is not set for the bs (fl == 0) and the strand is not minus.
        j1_tmp = j1;
        fl = 1;
      } else if (!fl && !strcmp(tag, "tes") && bs->strand != '+') { //if j1_tmp is not set for the bs (fl == 0) and the strand is not minus.
        j1_tmp = j1;
        fl = 1;
      }
      
      val_tmp = 0.0;
      for (;j1 ; j1 = j1->next) {
        if (j1->st >= ed) break; //if the sig pos is out of the win
        if (st > j1->st) tmp_st = st; //if st of sig block is up-stream pos of st
        else tmp_st = j1->st;
        if (j1->ed > ed) tmp_ed = ed; //if ed of sig block is down-stream pos of ed
        else tmp_ed = j1->ed;
        val_tmp += (j1->val) * (tmp_ed - tmp_st); //adding the val*len of sig block
      }

      if (chr_block_headsig_d) { //signal for denominator
        if (j1_tmp_d == NULL) j1_d = ch_sig_d->sig_list;
        else j1_d = j1_tmp_d;

        for (; j1_d; j1_d = j1_d->next) { //here's the slowest part...
          if (st < j1_d->ed && j1_d->st < ed) {
            break; //if one of sig block is inside the win
          } else if (j1_d->st >= ed) { //if there's no chance for j1 to overlap win
            j1_d = NULL;
            break;
          }
        }

        if (j1_d == NULL) { //if the win is the right side of the most right sig block
          printf("warning: signal denominator for region %lu-%lu on %s is zero. NA is returned.\n", st, ed, ch_smt->chr);
          if (add_one_val(ga_line_out, bs->line, "NA\n") != 0){
            LOG("error: output line was too long.");
            return -1; //making output link list with NA
          }
          ga_output_append (output_head, ga_line_out);
          continue;
        }

        if (!fl_d && !strcmp(tag, "tss") && bs->strand != '-') { //if j1_tmp is not set for the bs (fl_d == 0) and the strand is not minus.
          j1_tmp_d = j1_d;
          fl_d = 1;
        } else if (!fl_d && !strcmp(tag, "tes") && bs->strand != '+') { //if j1_tmp is not set for the bs (fl_d == 0) and the strand is not minus.
          j1_tmp_d = j1_d;
          fl_d = 1;
        }
        
        val_tmp_d = 0.0;
        for (;j1_d ; j1_d = j1_d->next) {
          if (j1_d->st >= ed) break; //if the sig pos is out of the win
          if (st > j1_d->st) tmp_st = st; //if st of sig block is up-stream pos of st
          else tmp_st = j1_d->st;
          if (j1_d->ed > ed) tmp_ed = ed; //if ed of sig block is down-stream pos of ed
          else tmp_ed = j1_d->ed;
          val_tmp_d += (j1_d->val) * (tmp_ed - tmp_st); //adding the val*len of sig block
        }
      }

      if (chr_block_headsig_d) {
        if (!val_tmp_d) {
          printf("warning: signal denominator for region %lu-%lu on %s is zero. NA is returned.\n", st, ed, ch_smt->chr);
          sprintf(tmp, "NA\n");
        } else {
//          sprintf(tmp, "%f\n", (val_tmp / val_tmp_d) / (float)(ed - st));
          sprintf(tmp, "%f\n", (val_tmp / val_tmp_d) );
        }
      }
      else sprintf(tmp, "%f\n", val_tmp / (float)(ed - st));

      if (add_one_val(ga_line_out, bs->line, tmp) != 0){
        LOG("error: output line was too long.");
        return -1; //making output link list with NA
      }
      ga_output_append (output_head, ga_line_out);
    } //bs
  } //chr

  return 0;
}
