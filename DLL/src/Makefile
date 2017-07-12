# Makefile for SRATS2010
#
#  Probably Rtools on cygwin is required to make it.
#

TARGET  = SRATS2010
DEF     = srats2010.def
DLLOBJS    = emsrm_windll.o
EMOBJS  = numlib.o em_exp.o em_gamma.o em_pareto.o em_tnorm.o em_lnorm.o em_tlogist.o em_llogist.o em_txvmax.o em_lxvmax.o em_txvmin.o em_lxvmin.o
OBJS = $(DLLOBJS)
OBJS := $(EMOBJS:%=emsrm/src/%)

CC      = gcc
DLLWRAP = dllwrap
CFLAGS  = -Wall -O2 -I emsrm/include
LDFLAGS = 

.SUFFIXES: .c .o

all: $(TARGET)

$(TARGET): $(OBJS) $(DEF)
	$(DLLWRAP) $(LDFLAGS) --dllname $@.dll --def $(DEF) $(DLLOBJS) $(EMOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(DLLOBJS) $(EMOBJS) $(TARGET).dll

