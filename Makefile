# Makefile for Canalysis
#

CC = gcc

#UNAME := $(shell uname -s)
#ifeq ($(UNAME),Darwin)
#	CC = gcc-mp-4.2
#endif

HDRS = util.h stats.h memalloc.h
UOBJS = libutil.o memalloc.o
OBJS = $(UOBJS) libstats.o

CFLAGS =  -g -Wall -O3 -std=c99
CPPFLAGS = -I$(HOME)/usr/include -I/opt/local/include
LDFLAGS += -lm -lcfitsio -lgsl -lgslcblas

INSTALL = /usr/bin/install
BINDIR = ${HOME}/usr/bin
TEST = ./test.pl

BINARIES = erf nsigma bin cstats rmsSlide smoothC lombscargle \
	pdm notchFilter getPeaks pca normC fftSmooth \
	fitsplitN2 fitsTSdump fitsTSsample ran_poisson \
	convolve text2fits fits2text bin2d interp text2fitsTS fasthdr

all: $(OBJS) $(BINARIES)

erf: erf.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@
nsigma: nsigma.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@

fasthdr: fasthdr.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) -o $@

# util only 
smoothC: smoothC.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

lombscargle: lombscargle.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

getPeaks: getPeaks.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

pca: pca.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

normC: normC.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

fftSmooth: fftSmooth.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

fitsplitN2: fitsplitN2.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

fitsTSdump: fitsTSdump.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

fitsTSsample: fitsTSsample.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

ran_poisson: ran_poisson.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

convolve: convolve.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

text2fits: text2fits.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

text2fitsTS: text2fitsTS.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

fits2text: fits2text.o $(UOBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(UOBJS) $(LDFLAGS) -o $@

interp: interp.o $(UOBJS) interpolate.o $(HDRS) interpolate.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< interpolate.o $(UOBJS) $(LDFLAGS) -o $@



# util + stats
bin: bin.o $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

bin2d: bin2d.o $(OBJS)  $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

cstats: cstats.o $(OBJS)  $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

pdm: pdm.o $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

notchFilter: notchFilter.o $(OBJS)  $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@

rmsSlide: rmsSlide.o $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(OBJS) $(LDFLAGS) -o $@


# local libs
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<


.PHONY: test clean install uninstall

test:
	$(TEST) $(BINARIES)

clean: 
	$(RM) $(OBJS) $(BINARIES) *.o core.* core *~ gmon.out *.testout

install:
	$(INSTALL) $(BINARIES) $(BINDIR)

uninstall:
	cd $(BINDIR); $(RM) $(BINARIES); cd -
