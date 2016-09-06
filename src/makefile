IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

EXEDIR = ../

_DEPS = ccs.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = matrixsize.o readgene.o bicluster.o bicluster_pair_score.o merge_bicluster.o print_bicluster.o ccs.o  
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEDIR)ccsbc: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
