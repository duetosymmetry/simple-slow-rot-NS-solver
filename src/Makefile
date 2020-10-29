IDIR = .
CC=clang
CXX=clang++
LD=clang++
BASEFLAGS=-I$(IDIR) $(shell gsl-config --cflags) -O3 \
  -Wall -Wno-narrowing
# if debugging, comment out this line:
BASEFLAGS += -DNOASSERT

CFLAGS   = $(BASEFLAGS) -std=c11
CXXFLAGS = $(BASEFLAGS) -std=c++1z

LDFLAGS  = -lm $(shell gsl-config --libs)

# this should contain all C and C++ source files
SRCS = ppEOS.cpp ppEOSTable.cpp util.cpp \
  testppEOS_eps_of_P.cpp testppEOS_get_rho_low_high.cpp \
  testppEOS.cpp testppEOS-cmdline.c \
  BackgroundModel.cpp \
  backgroundModelTest.cpp backgroundModelTest-cmdline.c \
  ExpAlpha0Beta0.cpp \
  testModelO1-single.cpp testModelO1-single-cmdline.c \
  testModelO1-multi.cpp testModelO1-multi-cmdline.c \
  writeModels.cpp ModelO2.cpp \
  runModelO2-multi-cmdline.c runModelO2-multi.cpp

TARGETS = testppEOS testppEOS_eps_of_P testppEOS_get_rho_low_high \
  backgroundModelTest \
  runPcPhicGrid \
	testPhi0Solver \
	runPhi0Family \
  testModelO1-single testModelO1-multi \
  runModelO2-multi \
  $(NULL)

.PHONY: all clean depend

# executable targets
all: $(TARGETS)

%-cmdline.o: %-cmdline.c %-cmdline.h

%-cmdline.c %-cmdline.h: %.ggo
	gengetopt < $<

# How to link/compile
$(TARGETS):
	$(LD) -o $@ $^ $(LDFLAGS)

# dependancies
testppEOS: testppEOS.o ppEOS.o ppEOSTable.o testppEOS-cmdline.o

testppEOS_eps_of_P: testppEOS_eps_of_P.o ppEOS.o

testppEOS_get_rho_low_high: testppEOS_get_rho_low_high.o ppEOS.o

backgroundModelTest: backgroundModelTest.o                           \
                     backgroundModelTest-cmdline.o BackgroundModel.o \
										 ExpAlpha0Beta0.o \
                     writeModels.o ppEOS.o ppEOSTable.o

runPcPhicGrid: runPcPhicGrid.o                           \
               runPcPhicGrid-cmdline.o BackgroundModel.o \
						   ExpAlpha0Beta0.o \
							 util.o \
               writeModels.o ppEOS.o ppEOSTable.o

testPhi0Solver: testPhi0Solver.o \
                testPhi0Solver-cmdline.o \
                BackgroundModel.o \
								ExpAlpha0Beta0.o \
                Phi0Solver.o \
                writeModels.o ppEOS.o ppEOSTable.o

runPhi0Family: runPhi0Family.o \
               runPhi0Family-cmdline.o \
               BackgroundModel.o \
							 ExpAlpha0Beta0.o \
               Phi0Solver.o \
               writeModels.o ppEOS.o ppEOSTable.o

testModelO1-single: testModelO1-single.o testModelO1-single-cmdline.o\
                    BackgroundModel.o ModelO1.o writeModels.o        \
                    ppEOS.o ppEOSTable.o ModelO2.o

testModelO1-multi: testModelO1-multi.o testModelO1-multi-cmdline.o \
                   BackgroundModel.o ModelO1.o writeModels.o       \
                   ppEOS.o ppEOSTable.o ModelO2.o

runModelO2-multi: runModelO2-multi.o runModelO2-multi-cmdline.o     \
                  BackgroundModel.o ModelO1.o writeModels.o ppEOS.o \
                  ppEOSTable.o ModelO2.o util.o

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
	makedepend $(DEPFLAGS) -- $(BASEFLAGS) -- $(SRCS) >&/dev/null

# put this file in the last line of your Makefile
sinclude $(DEPFILE)
# DO NOT DELETE
