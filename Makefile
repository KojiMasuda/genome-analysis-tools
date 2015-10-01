CC=gcc
TARGET=ga_overlap
OBJS=ga_overlap.o parse_chr.o write_tab.o argument.o
CFLAGS+=-O
CFLAGS+=-Wall
.SUFFIXES: .c .o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(CFLAGS)

.c.o:  $<
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f $(OBJS) $(TARGET)
