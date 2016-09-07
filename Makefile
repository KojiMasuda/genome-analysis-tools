CC=gcc
OBJS1=ga_overlap.o parse_chr.o write_tab.o argument.o
OBJS2=ga_reads_summit.o parse_chr.o write_tab.o argument.o sort_list.o ga_math.o
OBJS3=ga_reads_summit_all.o parse_chr.o write_tab.o argument.o sort_list.o ga_math.o
OBJS4=ga_calc_dist.o parse_chr.o write_tab.o argument.o sort_list.o
OBJS5=ga_reads_region.o parse_chr.o write_tab.o argument.o sort_list.o
OBJS6=ga_deltaG.o parse_chr.o write_tab.o argument.o
OBJS7=ga_nuc_region.o parse_chr.o write_tab.o argument.o sort_list.o

TARGET=ga_overlap ga_reads_summit ga_reads_summit_all ga_calc_dist ga_reads_region ga_deltaG ga_nuc_region
#CFLAGS+=-O3
CFLAGS+=-O0
CFLAGS+=-g
CFLAGS+=-Wall
LIBS += -lz -lm
.SUFFIXES: .c .o

all: ga_overlap ga_reads_summit ga_reads_summit_all ga_calc_dist ga_reads_region ga_deltaG ga_nuc_region

ga_overlap: $(OBJS1)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_reads_summit: $(OBJS2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_reads_summit_all: $(OBJS3)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_calc_dist: $(OBJS4)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_reads_region: $(OBJS5)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_deltaG: $(OBJS6)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_nuc_region: $(OBJS7)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.c.o:  $<
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6) $(OBJS7) $(TARGET)
