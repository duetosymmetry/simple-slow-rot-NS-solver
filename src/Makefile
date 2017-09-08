IDIR = .
CXX=clang++
CC=$(CXX)
LD=clang++
CFLAGS=-I$(IDIR) $(shell gsl-config --cflags) -O3 \
  -Wall -Wno-narrowing -std=c++1z
# if debugging, comment out this line:
CFLAGS += -DNOASSERT
CXXFLAGS=$(CFLAGS)

LDFLAGS=-lm $(shell gsl-config --libs)

_DEPS = parameters.h declaration.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# this should contain all C and C++ source files
SRCS = ppEOS.cc testppEOS.cc testppEOS_eps_of_P.cc                    \
  testppEOS_get_rho_low_high.cc ppEOSTable.cc BackgroundModel.cc      \
  backgroundModelTest-cmdline.c testppEOS-cmdline.c                   \
  testModelO1-single-cmdline.c testModelO1-multi-cmdline.c ModelO1.cc \
  testModelO1-single.cc testModelO1-multi.cc ModelO2.cc               \
  runModelO2-multi-cmdline.c runModelO2-multi.cc util.cc

TARGETS = testppEOS testppEOS_eps_of_P testppEOS_get_rho_low_high \
  backgroundModelTest testModelO1-single testModelO1-multi        \
  runModelO2-multi

.PHONY: all clean depend

# executable targets
all: $(TARGETS)

%-cmdline.o: %-cmdline.c %-cmdline.h

%-cmdline.c %-cmdline.h: %.ggo
	gengetopt < $<

testppEOS: testppEOS.o ppEOS.o ppEOSTable.o testppEOS-cmdline.o
testppEOS_eps_of_P: testppEOS_eps_of_P.o ppEOS.o
testppEOS_get_rho_low_high: testppEOS_get_rho_low_high.o ppEOS.o

backgroundModelTest: backgroundModelTest.o                           \
                     backgroundModelTest-cmdline.o BackgroundModel.o \
                     writeModels.o ppEOS.o ppEOSTable.o              \
                     ModelO1.o ModelO2.o
	$(LD) -o $@ $^ $(LDFLAGS)

testModelO1-single: testModelO1-single.o testModelO1-single-cmdline.o\
                    BackgroundModel.o ModelO1.o writeModels.o        \
                    ppEOS.o ppEOSTable.o ModelO2.o
	$(LD) -o $@ $^ $(LDFLAGS)

testModelO1-multi: testModelO1-multi.o testModelO1-multi-cmdline.o \
                   BackgroundModel.o ModelO1.o writeModels.o       \
                   ppEOS.o ppEOSTable.o ModelO2.o
	$(LD) -o $@ $^ $(LDFLAGS)

runModelO2-multi: runModelO2-multi.o runModelO2-multi-cmdline.o     \
                  BackgroundModel.o ModelO1.o writeModels.o ppEOS.o \
                  ppEOSTable.o ModelO2.o util.o
	$(LD) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *~ core $(INCDIR)/*~ $(TARGETS)


###
# makedepends related settings
# -Y    Don't search standard library paths, this is excessive
# -f    output dependencies to this file
# -s    use this token string
###

DEPFILE   = .depends
DEPTOKEN  = '\# MAKEDEPENDS'
DEPFLAGS  = -Y -f $(DEPFILE) -s $(DEPTOKEN)

depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

$(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRCS) >&/dev/null

# put this file in the last line of your Makefile
sinclude $(DEPFILE)
# DO NOT DELETE
