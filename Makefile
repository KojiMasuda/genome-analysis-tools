CC=gcc
OBJS1=ga_overlap.o parse_chr.o write_tab.o argument.o
OBJS2=ga_reads_summit.o parse_chr.o write_tab.o argument.o sort_list.o ga_math.o

TARGET=ga_overlap ga_reads_summit
#CFLAGS+=-O3
CFLAGS+=-O0
CFLAGS+=-g
CFLAGS+=-Wall
LIBS += -lz -lm
.SUFFIXES: .c .o

all: ga_overlap ga_reads_summit

ga_overlap: $(OBJS1)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
ga_reads_summit: $(OBJS2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.c.o:  $<
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET)
