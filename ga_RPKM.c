/*
 * This program is one of the genome analysis tools.
 * This program calculates expression level as reads per kilobase per million mapped reads (RPKM). Input file is normalized mapped reads (per million) in bedgraph format.
 * See usage for detail.
 */

#include "parse_chr.h"
#include "write_tab.h"
#include "sort_list.h"
#include "argument.h"
#include "ga_my.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define LOG(m) \
  fprintf(stderr, \
  "%s:line%d:%s(): " m "\n", \
  __FILE__, __LINE__, __FUNCTION__)

static int make_exon_arr (unsigned long arr[], char *exon);
static int reassign_exon (struct ref *ref);

static void usage()
{
  printf("Tool:    ga_RPKM\n\n\
Summary: report the expression level as reads per kilobase per million mapped reads (RPKM)\n\n\
Usage:   ga_RPKM [options] -exp <expression file in bedgraph format> -ref <reference file> -readlen <read length(default:101)>\n\
or       ga_RPKM [options] -exp <expression file in bedgraph format> -ref <reference file> -readlen <read length(default:101)> --consid_ov --thresh <threshold(default:0.5)>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --col_chr <int>: column number for chromosome of reference file (default:2).\n\
         --col_start <int>: column number for gene start position of reference file (default:4).\n\
         --col_end <int>: column number for gene end position of reference file (default:5).\n\
         --col_strand <int>: column number for strand of reference file (default:3).\n\
         --col_exon_st <int>: column number for exon start positions of reference file (default:9).\n\
         --col_exon_ed <int>: column number for exon end positions of reference file (default:10).\n\
         --col_gene <int>: column number for gene name of reference file (default:0).\n\
         --consid_ov: consider exon overlapping. If exon is overlapped less than (1-thresh), the non-overlapping exon is used for calculation. Set --thresh argument(default: off)l\n\
         --thresh: threshold for considering exon overlapping. The 'non-overlapping' exon is used for calculation for exon which is overlapped by proportion of '1-thresh'(default:0.5).\n\
         --header: the first line of reference file is considered as header (default:off).\n");
  exit(0);
}

static void version()
{
  printf("'ga_RPKM' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}


static int hf = 0;
static char hfs[4] = "off\0";
static int cf = 0;
static char cfs[4] = "off\0";
static char *fileexp = NULL;
static char *fileref = NULL;
static int readlen = 101;
static int col_chr = 2;
static int col_st = 4;
static int col_ed = 5;
static int col_strand = 3;
static int col_ex_st = 9;
static int col_ex_ed = 10;
static int col_gene = 0;
static double thresh = 0.5; //threshold for considering exon overlapping. The "non-overlapping" exon is used for calculation for exon of which is overlapped proportion is less than (1-thresh). Thus if thresh is 0.7, and overlapping rate is less than 0.3, the remaining non-overlapping (>=0.7) exon is used for calculation.
char *ga_header_line = NULL; //header line. Note this is external global variable
static char ga_line_out[LINE_STR_LEN] = {0}; //output line

static const Argument args[] = {
  {"-h"           , ARGUMENT_TYPE_FUNCTION, usage        },
  {"--help"       , ARGUMENT_TYPE_FUNCTION, usage        },
  {"-v"           , ARGUMENT_TYPE_FUNCTION, version      },
  {"--header"     , ARGUMENT_TYPE_FLAG_ON , &hf          },
  {"--consid_ov"  , ARGUMENT_TYPE_FLAG_ON , &cf          },
  {"-exp"         , ARGUMENT_TYPE_STRING  , &fileexp     },
  {"-ref"         , ARGUMENT_TYPE_STRING  , &fileref     },
  {"-readlen"     , ARGUMENT_TYPE_INTEGER , &readlen     },
  {"--col_chr"    , ARGUMENT_TYPE_INTEGER , &col_chr     },
  {"--col_start"  , ARGUMENT_TYPE_INTEGER , &col_st      },
  {"--col_end"    , ARGUMENT_TYPE_INTEGER , &col_ed      },
  {"--col_strand" , ARGUMENT_TYPE_INTEGER , &col_strand  },
  {"--col_exon_st", ARGUMENT_TYPE_INTEGER , &col_ex_st   },
  {"--col_exon_ed", ARGUMENT_TYPE_INTEGER , &col_ex_ed   },
  {"--col_gene"   , ARGUMENT_TYPE_INTEGER , &col_gene    },
  {"--thresh"     , ARGUMENT_TYPE_FLOAT   , &thresh      },
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
  unsigned long arr_ex_st[1000]; //array of exon positions
  unsigned long arr_ex_ed[1000];
  unsigned long st;
  unsigned long tmp_st;
  unsigned long ed;
  unsigned long tmp_ed;
  unsigned long exon_len;
  int fl, c, exon_nb;
  double val;

  time_t timer;

  if(hf) strcpy(hfs, "on\0");
  if(cf) strcpy(cfs, "on\0");
  time(&timer);
  printf("Tool:                            %s\n\n\
Input expression file:             %s\n\
Input reference file:              %s\n\
Read length:                       %d\n\
ref col of chr, start, end, strand:%d, %d, %d, %d\n\
ref col exon start, end:           %d, %d\n\
ref col gene:                      %d\n\
header flag:                       %s\n\
consider overlapping?:             %s\n\
threshold:                         %.3f\n\
time:                              %s\n",\
 "ga_RPKM", fileexp, fileref, readlen, col_chr, col_st, col_ed, col_strand, col_ex_st, col_ex_ed, col_gene, hfs, cfs, thresh, ctime(&timer) );

  ga_parse_file_path (fileexp, path_exp, fn_exp, ext_exp); //parsing input file name into path, file name, and extension
  ga_parse_file_path (fileref, path_ref, fn_ref, ext_ref);

  ga_parse_chr_ref (fileref, &chr_block_headref, col_chr, col_st, col_ed, col_strand, col_ex_st, col_ex_ed, col_gene, hf); //parsing each ref for each chromosome

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

    if (cf) { //if considering exon overlap
      if(reassign_exon (ch_ref -> ref_list) != 0) {
        LOG("error: error in reassign_exon function.");
        goto err;
      }
    }

    if (ch_exp == NULL) {
      for (r = ch_ref -> ref_list; r; r = r -> next) {
        if (cf) sprintf(tmp, "0\t%s\t%s\t%s\t%.3f\n", r -> rm_ex_st, r -> rm_ex_ed, r -> ov_gene, r -> ov_prop);
        else sprintf(tmp, "0\n");

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
      val = 0;
      exon_len = 0;

      if (cf) {
        if (0 < r -> ov_prop && r -> ov_prop < 1 - thresh) { //if partial overlapping
          exon_nb = make_exon_arr (arr_ex_st, r -> rm_ex_st);
          exon_nb = make_exon_arr (arr_ex_ed, r -> rm_ex_ed);
        } else { //if no overlapping or overlapping rate is more than 1 - thresh
          exon_nb = make_exon_arr (arr_ex_st, r -> ex_st);
          exon_nb = make_exon_arr (arr_ex_ed, r -> ex_ed);
        }
      } else {
        exon_nb = make_exon_arr (arr_ex_st, r -> ex_st);
        exon_nb = make_exon_arr (arr_ex_ed, r -> ex_ed);
      }

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

        if (!(r -> ov_prop)) { //if no overlap
          if (j1_tmp == NULL) { //if j1_tmp is not set for the chr
            j1_tmp = j1;
            fl = 1;
          } else if (!fl) { //else if j1_tmp is not set for the ref (fl == 0)
            j1_tmp = j1;
            fl = 1;
          }
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

      if (cf) sprintf(tmp, "%.6f\t%s\t%s\t%s\t%.3f\n", val / ( (float)exon_len / 1000 ) / (float)readlen, r -> rm_ex_st, r -> rm_ex_ed, r -> ov_gene, r -> ov_prop);
      else sprintf(tmp, "%.6f\n", val / ( (float)exon_len / 1000 ) / (float)readlen ); //calculating RPKM
      if (add_one_val(ga_line_out, r -> line, tmp) < 0) {
        LOG("error: output line was too long.");
        goto err;
      }
      ga_output_append (&output_head, ga_line_out);
    } //for ref
  } //for chr

  if (cf) sprintf(output_name, "%sRPKM_of_%s_%s_consid_ov.txt", path_exp, fn_exp, fn_ref);
  else sprintf(output_name, "%sRPKM_of_%s_%s.txt", path_exp, fn_exp, fn_ref);

  if (hf && !cf) add_one_val(ga_line_out, ga_header_line, "RPKM\n");
  else if (hf && cf) add_one_val(ga_line_out, ga_header_line, "RPKM\tsub_exons_st\tsub_exons_ed\tov_genes\tov_prop\n");
  if (hf) ga_write_lines (output_name, output_head, ga_line_out);
  else ga_write_lines (output_name, output_head, ga_header_line);

  if (chr_block_headref) ga_free_chr_block(&chr_block_headref);
  if (chr_block_headexp) ga_free_chr_block(&chr_block_headexp);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) my_free(ga_header_line);
  return 0;

err:
  if (chr_block_headref) ga_free_chr_block(&chr_block_headref);
  if (chr_block_headexp) ga_free_chr_block(&chr_block_headexp);
  if (output_head) ga_free_output(&output_head);
  if (ga_header_line) my_free(ga_header_line);
  return -1;
}

/*
 * this makes array of exon positions and returns exon number
 * arr: array of exon positions
 * exon: char pointer to exon such as "123,456,".
 */
static int make_exon_arr (unsigned long arr[], char *exon)
{
  char ex[80] = {0}; //one exon position like "123".
  char *end;
  int c, e, e_tmp;

  e = 0; //char position for exon
  e_tmp = 0; //char position for tmp exon
  c = 0; //number of exon
  while(exon[e_tmp] != '\0') {
    if (exon[e_tmp] == ',') {
      ex[e] = '\0';
      arr[c] = strtoul(ex, &end, 0);
      e = 0;
      e_tmp++;
      c++;
      continue;
    } else {
      ex[e] = exon[e_tmp]; //extracting start pos of each exon
    }
    e++;
    e_tmp++;
  } //exon start

  return c;
}

static int reassign_exon (struct ref *ref)
{

  struct ref *i, *j;
  unsigned long arr_exi_st[1000];
  unsigned long arr_exi_ed[1000];
  unsigned long arr_exj_st[1000];
  unsigned long arr_exj_ed[1000];
  unsigned long exon_len, rm_exon_len; //exon length and remaining exon length
  unsigned long ed_tmp;
  int ci, cj, exoni_nb, exonj_nb, gene_ov, exon_ov, co;
  char *tmp;
  char ex_st_tmp[80] = {0}, ex_ed_tmp[80] = {0};

  for (i = ref; i; i = i -> next) {
    exoni_nb = make_exon_arr(arr_exi_st, i -> ex_st);
    exoni_nb = make_exon_arr(arr_exi_ed, i -> ex_ed);
    gene_ov = 0; //gene overlapping flag
    exon_ov = 0; //exon overlapping flag
    co      = 0; //complete overlapping flag
    exon_len = 0; //init exon len
    rm_exon_len = 0; //init remaining exon len

    for (ci = 0; ci < exoni_nb; ci++) { //calculating original exon length
      exon_len += arr_exi_ed[ci] - arr_exi_st[ci] + 1;
    }

    for (j = ref; j; j = j -> next) {
      if (i == j) continue; //if i == j
      if (i->st < j->ed && j->st < i->ed) { //if genes are overlapping
        gene_ov = 1; //flag is on
        if (NULL == i->ov_gene) { //if the overlapping is the first one
          i -> ov_gene = (char *)my_malloc(sizeof(char) * (strlen(j -> gene) + 2)); //"gene" + ',' + '\0'
          sprintf(i -> ov_gene, "%s,", j -> gene);
        } else { //if the overlapping is not the first one
          tmp = (char *)my_realloc(i -> ov_gene, sizeof(char) * (strlen(i -> ov_gene) + strlen(j -> gene) + 2)); //"gene" + ',' + '\0'
          i -> ov_gene = tmp;
          sprintf(i -> ov_gene, "%s%s,", i -> ov_gene, j -> gene);
        }

        exonj_nb = make_exon_arr(arr_exj_st, j -> ex_st);
        exonj_nb = make_exon_arr(arr_exj_ed, j -> ex_ed);
        for (ci = 0; ci < exoni_nb; ci++) {
          for (cj = 0; cj < exonj_nb; cj++) {
            if (arr_exi_st[ci] < arr_exj_ed[cj] && arr_exj_st[cj] < arr_exi_ed[ci] ) { //if exons are overlapping
              exon_ov = 1; //flag is on
              if (arr_exi_st[ci] > arr_exj_st[cj] && arr_exi_ed[ci] > arr_exj_ed[cj] ) { //if i exon partially overlaps j exon
                arr_exi_st[ci] = arr_exj_ed[cj];
              } else if (arr_exi_st[ci] >= arr_exj_st[cj] && arr_exi_ed[ci] <= arr_exj_ed[cj] ) { //if i exon is completely overlapped by j exon
                if(sizeof(arr_exi_st)/sizeof(arr_exi_st[0]) < exoni_nb) {
                  LOG("error: exon number exceeds the exon array.");
                  goto err;
                }
                arr_delete_ul(arr_exi_st, exoni_nb, ci);

                if(sizeof(arr_exi_ed)/sizeof(arr_exi_ed[0]) < exoni_nb) {
                  LOG("error: exon number exceeds the exon array.");
                  goto err;
                }
                arr_delete_ul(arr_exi_ed, exoni_nb, ci);

                co = 1; //complete overlapping flag on
                exoni_nb--;
              } else if (arr_exi_st[ci] < arr_exj_st[cj] && arr_exi_ed[ci] < arr_exj_ed[cj] ) { //if i exon partially overlaps j exon
                arr_exi_ed[ci] = arr_exj_st[cj];
              } else if (arr_exi_st[ci] < arr_exj_st[cj] && arr_exi_ed[ci] > arr_exj_ed[cj] ) { //if i exon completely covers j exon
                ed_tmp = arr_exi_ed[ci]; //preserving ed position
                arr_exi_ed[ci] = arr_exj_st[cj];

                if(sizeof(arr_exi_st)/sizeof(arr_exi_st[0]) < exoni_nb + 1) {
                  LOG("error: exon number exceeds the exon array.");
                  goto err;
                }
                arr_insert_ul(arr_exi_st, exoni_nb, ci+1, arr_exj_ed[cj]); //inserting new exon

                if(sizeof(arr_exi_ed)/sizeof(arr_exi_ed[0]) < exoni_nb + 1) {
                  LOG("error: exon number exceeds the exon array.");
                  goto err;
                }
                arr_insert_ul(arr_exi_ed, exoni_nb, ci+1, ed_tmp); //inserting new exon

                exoni_nb++;
              }
            }//if exon overlapping
          }//for each j exon
          if (co) {
            ci--;
            co = 0;
          }
        }//for each i exon
      } //if gene overlapping
    } //for j(gene)
    if (!gene_ov) i -> ov_gene = strdup("NA");

    if(exon_ov) {
      if (!exoni_nb) { //if all exons are overlapped...
        i -> rm_ex_st = strdup("NA");
        i -> rm_ex_ed = strdup("NA");
        i -> ov_prop = 1;
      } else {
        for (ci = 0; ci < exoni_nb; ci++) { //appending remaining exons and calculating remaining exon length
          sprintf(ex_st_tmp, "%lu,", arr_exi_st[ci]);
          sprintf(ex_ed_tmp, "%lu,", arr_exi_ed[ci]);

          if (NULL == i->rm_ex_st) { //if the exon is the first one
            i -> rm_ex_st = (char *)my_malloc(sizeof(char) * (strlen(ex_st_tmp) + 1)); //"exon_st" + ',' + '\0'
            i -> rm_ex_ed = (char *)my_malloc(sizeof(char) * (strlen(ex_ed_tmp) + 1)); //"exon_ed" + ',' + '\0'
            sprintf(i -> rm_ex_st, "%s", ex_st_tmp);
            sprintf(i -> rm_ex_ed, "%s", ex_ed_tmp);
          } else { //if the exon is not the first one
            tmp = (char *)my_realloc(i -> rm_ex_st, sizeof(char) * (strlen(i -> rm_ex_st) + strlen(ex_st_tmp) + 1)); //"exon_st" + ',' + '\0'
            i -> rm_ex_st = tmp;
            sprintf(i -> rm_ex_st, "%s%s", i -> rm_ex_st, ex_st_tmp);
            tmp = (char *)my_realloc(i -> rm_ex_ed, sizeof(char) * (strlen(i -> rm_ex_ed) + strlen(ex_ed_tmp) + 1)); //"exon_ed" + ',' + '\0'
            i -> rm_ex_ed = tmp;
            sprintf(i -> rm_ex_ed, "%s%s", i -> rm_ex_ed, ex_ed_tmp);
          }

          rm_exon_len += arr_exi_ed[ci] - arr_exi_st[ci] + 1;
        }
        i -> ov_prop = (exon_len - rm_exon_len) / (double)exon_len;
      }
    } else {
      i -> rm_ex_st = strdup(i -> ex_st);
      i -> rm_ex_ed = strdup(i -> ex_ed);
      i -> ov_prop = 0;
    }
  } //for i(gene)

  return 0;

err:

  return -1;
}

