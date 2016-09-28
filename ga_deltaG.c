/*
 * This program is one of the genome analysis tools.
 * This program makes the wiggle file of the free energy difference between the duplex and single-strand states.
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


static void usage()
{
  printf("Tool:    genome analysis\n\n\
Program: ga_deltaG\n\
Summary: makes the wiggle file of the free energy(delta-G) difference between the duplex and single-strand states.\n\n\
Usage:   ga_deltaG -fa <fasta file> -gt <genome table file> -lib <library name: SantaLucia_1996 | Breslauer_1986> -win <window size(int)> -step <step size(int)>\n\n\
Options:\n\
         -v: output version information and exit.\n\
         -h, --help: display this help and exit.\n\
         --gradient: coefficients are used for the energy calculation. The middle of the window is coeffient of 1, and it linearly decreases to 0 for both left and right side (default:0).\n");
  exit(0);
}

static void version()
{
  printf("'ga_deltaG' in genome analysis tools version: %d.%d.%d\n", VER_MAJOR, VER_MOD, VER_MINOR);
  exit(0);
}

char *ga_header_line = NULL; //header line. Note this is external global variable
static char *fa = NULL; //fasta
static char *gt = NULL; //genome table
static char *lib = NULL; //library
static int step = 0, win = 0; //step and window size
static int gf = 0; //gradient flag


static const Argument args[] = {
  {"-h"          , ARGUMENT_TYPE_FUNCTION, usage  },
  {"--help"      , ARGUMENT_TYPE_FUNCTION, usage  },
  {"-v"          , ARGUMENT_TYPE_FUNCTION, version},
  {"-fa"         , ARGUMENT_TYPE_STRING  , &fa    },
  {"-gt"         , ARGUMENT_TYPE_STRING  , &gt    },
  {"-lib"        , ARGUMENT_TYPE_STRING  , &lib   },
  {"-step"       , ARGUMENT_TYPE_INTEGER , &step  },
  {"-win"        , ARGUMENT_TYPE_INTEGER , &win   },
  {"--gradient"  , ARGUMENT_TYPE_FLAG_ON , &gf    },
  {NULL          , ARGUMENT_TYPE_NONE    , NULL   },
};

int main (int argc, char *argv[])
{
  /*path, filename, and extension*/
  char path[PATH_STR_LEN] = {0};
  char fn[FILE_STR_LEN] = {0};
  char ext[EXT_STR_LEN] = {0};
  char output_name[PATH_STR_LEN + FILE_STR_LEN + EXT_STR_LEN] = {0}; //output file name
  char tmp[LINE_STR_LEN] = {0};
  char *frag; //fragment
  char di[3] = {0}; //di-Nucleotide

  struct chr_block_fa *chr_block_head = NULL; //for parsing fasta
  struct chr_block_fa *ch; //for
  struct chr_block *chr_block_head_gt = NULL; //for parsing genome table

  float g_AA, g_AT, g_TA, g_CA, g_GT, g_CT, g_GA, g_CG, g_GC, g_GG, dG; //nearest-neighbor values and delta-G
  float co, grad; //coefficient and gradient
  unsigned long i;
  int s, mid;

  time_t timer;
  struct gzFile_s *gfp = NULL;

  argument_read(&argc, argv, args);//reading arguments
  if (fa == NULL || gt == NULL || lib == NULL || !win || !step) usage();

  time(&timer);
  printf("Program:                         %s\n\
Tools:                           genome analysis tools\n\n\
fasta file:                      %s\n\
genome table:                    %s\n\
library:                         %s\n\
window size:                     %d\n\
step size:                       %d\n\
gradient?:                       %d\n\
time:                            %s\n",\
 "ga_deltaG", fa, gt, lib, win, step, gf, ctime(&timer) );

  ga_parse_chr_bs (gt, &chr_block_head_gt, 0, 1, 1, -1, 0); //parsing genome tible
  if (ga_parse_chr_fa(fa, &chr_block_head, chr_block_head_gt) != 0){
    LOG("error: error in ga_parse_chr_fa function.");
    goto err;
  }
  ga_parse_file_path (fa, path, fn, ext); //parsing input file name into path, file name, and extension

  //for now, delta G library of SantaLucia_1996 and Breslauer_1986 are supported.
  if (!strcmp(lib, "SantaLucia_1996")) {
    g_AA = 1.02; //also for TT
    g_AT = 0.73;
    g_TA = 0.60;
    g_CA = 1.38; //also for TG
    g_GT = 1.43; //also for AC
    g_CT = 1.16; //also for AG
    g_GA = 1.46; //also for TC
    g_CG = 2.09;
    g_GC = 2.28;
    g_GG = 1.77; //also for CC
  } else if (!strcmp(lib, "Breslauer_1986")) {
    g_AA = 1.9; //also for TT
    g_AT = 1.5;
    g_TA = 0.9;
    g_CA = 1.9; //also for TG
    g_GT = 1.3; //also for AC
    g_CT = 1.6; //also for AG
    g_GA = 1.6; //also for TC
    g_CG = 3.6;
    g_GC = 3.1;
    g_GG = 3.1; //also for CC
  } else {
    printf ("Please specify correct delta G library: %s.\n", lib );
    goto err;
  }

  frag = (char *)malloc(sizeof(char) * (win + 1)); //fragment DNA from genome
  for (ch = chr_block_head; ch; ch = ch->next) {
    printf("calculating on %s \n", ch->chr);

    if (gf) sprintf(output_name, "%sdeltaG.%s.win%d_grad_lib-%s_%s.%d.wig.gz", path, fn, win, lib, ch->chr, step);
    else sprintf(output_name, "%sdeltaG.%s.win%d_lib-%s_%s.%d.wig.gz", path, fn, win, lib, ch->chr, step);
    if ((gfp = gzopen(output_name, "w")) == NULL) {
      LOG("error: output file cannot be open.");
      exit(EXIT_FAILURE);
    }

    sprintf(tmp, "variableStep\tchrom=%s\tspan=%d\n", ch->chr, step); //step size is used as "span" here...
    if (gzputs(gfp, tmp)<0){
      LOG("error: file writing error.");
      goto err;
    }

    if (gf) {
      mid = win/2 - 1; //the middle position of the window
      grad = 2.0 / (float)win; //the coefficient of gradient
    }
    for (i = 0; i < (ch->letter_len) / step - 1; i++) {
      if (ch->letter_len <= i*step + (win-1) ) break; //if window exceeds chromosome end...
      strncpy(frag, ch->letter + i*step, win);  //sliced fragment DNA from genome. Starting from letter[0], win nucleotide is copied to frag. Next the slice starts from letter[0 + i*step].
      frag[win] = '\0'; //null
      if (strchr(frag, 'N') != NULL) continue; //if N is found
      dG = 0.0; //init delta-G

      if (gf) { //if gradient flag is on
        co = 2.0 / (float)win; //the first coefficient

        for (s = 0; s < strlen(frag) - 1; s++) {
          strncpy(di, frag + s, 2); //extracting di-Nucleotide, di[2] is always '\0'.

          if (!strcmp(di, "AA") || !strcmp(di, "TT")) dG += g_AA * co; //also for TT
          else if (!strcmp(di, "AT")) dG += g_AT * co; //also for 
          else if (!strcmp(di, "TA")) dG += g_TA * co; //also for 
          else if (!strcmp(di, "CA") || !strcmp(di, "TG")) dG += g_CA * co; //also for 
          else if (!strcmp(di, "GT") || !strcmp(di, "AC")) dG += g_GT * co; //also for 
          else if (!strcmp(di, "CT") || !strcmp(di, "AG")) dG += g_CT * co; //also for 
          else if (!strcmp(di, "GA") || !strcmp(di, "TC")) dG += g_GA * co; //also for 
          else if (!strcmp(di, "CG")) dG += g_CG * co; //also for 
          else if (!strcmp(di, "GC")) dG += g_GC * co; //also for 
          else if (!strcmp(di, "GG") || !strcmp(di, "CC")) dG += g_GG * co; //also for 
          else printf ("Oops!! Your DNA sequence is not expected!!: %s.\n", di);

          if (s >= mid) co = co - grad; //decreasing co
          else co = co + grad; //increasing co until the mid position
        }
      } else {
        for (s = 0; s < strlen(frag) - 1; s++) {
          strncpy(di, frag + s, 2); //extracting di-Nucleotide, di[2] is always '\0'.

          if (!strcmp(di, "AA") || !strcmp(di, "TT")) dG += g_AA; //also for TT
          else if (!strcmp(di, "AT")) dG += g_AT; //also for 
          else if (!strcmp(di, "TA")) dG += g_TA; //also for 
          else if (!strcmp(di, "CA") || !strcmp(di, "TG")) dG += g_CA; //also for 
          else if (!strcmp(di, "GT") || !strcmp(di, "AC")) dG += g_GT; //also for 
          else if (!strcmp(di, "CT") || !strcmp(di, "AG")) dG += g_CT; //also for 
          else if (!strcmp(di, "GA") || !strcmp(di, "TC")) dG += g_GA; //also for 
          else if (!strcmp(di, "CG")) dG += g_CG; //also for 
          else if (!strcmp(di, "GC")) dG += g_GC; //also for 
          else if (!strcmp(di, "GG") || !strcmp(di, "CC")) dG += g_GG; //also for 
          else printf ("Oops!! Your DNA sequence is not expected!!: %s.\n", di);
        }
      }

      sprintf (tmp, "%lu\t%.2f\n", (i*step+(win/2)-(step/2)), dG);
      if (gzputs(gfp, tmp)<0){
        LOG("error: file writing error.");
        goto err;
      }
    } //each step
    gzclose(gfp);
  } //chromosome

  ga_free_chr_block_fa(&chr_block_head);
  ga_free_chr_block(&chr_block_head_gt);
  free(frag);
  return 0;

err:
  gzclose(gfp);
  if (chr_block_head) ga_free_chr_block_fa(&chr_block_head);
  if (chr_block_head_gt) ga_free_chr_block(&chr_block_head_gt);
  if (frag) free(frag);
  return -1;
}

