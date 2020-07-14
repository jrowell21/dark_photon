CXX=g++
CXXFLAGS=-O3 -fPIC -Wall

ROOT_FLAGS=$(shell $(ROOTSYS)/bin/root-config --cflags)
ROOT_INCLUDE=$(shell $(ROOTSYS)/bin/root-config  --incdir)

BUILDDIR=.
OBJ_DIR=.


ROOT_LIBS    = $(shell $(ROOTSYS)/bin/root-config --libs) -lTreePlayer -lMinuit -lXMLIO -lMLP -lRIO -lTMVA 
GLIBS := $(shell root-config --glibs)
FLAGS = $(CXXFLAGS)
FLAGS += ${ROOT_FLAGS}
FLAGS_ROOFIT = $(FLAGS)
FLAGS_ROOFIT += -lRooFitCore -lRooFit -lRooStats -lFoam

ROODCB=DSCB
ROODCBIR=../RooFit-pdfs


MACRO=mass_calibration

all: $(MACRO)

$(ROODCBIR)/$(ROODCB).o:$(ROODCBIR)/src/$(ROODCB).cc $(ROODCBIR)/include/$(ROODCB).h
	@echo "---> Making $(ROODCB)"
	$(CXX) $(CXXFLAGS) $(ROOT_LIBS) $(GLIBS) $(FLAGS_ROOFIT) -I $(ROODCBIR)/include -I $(ROOT_INCLUDE) -c $< -o $@

$(MACRO): $(ROODCBIR)/$(ROODCB).o  main.cpp 
	@echo "---> Making mass_calibration.exe..."
	$(CXX) $(CXXFLAGS) $(ROOT_LIBS) $(GLIBS) $(FLAGS_ROOFIT) -I $(ROODCBIR)/include -I $(ROOT_INCLUDE) $^ -o $(MACRO).exe


clean:
	rm -f $(BUILDDIR)/$(MACRO).exe $(ROODCBIR)/$(ROODCB).o
