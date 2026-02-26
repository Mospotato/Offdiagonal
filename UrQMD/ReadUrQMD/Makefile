ROOTFLAG = `root-config --cflags`
ROOTLIBS = `root-config --libs`

CODE=Code
DEPDIR = ../Dependency
VPATH = $(CODE):$(DEPDIR)
SRCS = $(shell find $(CODE) $(DEPDIR) -name "*.cpp" ! -name "KochRatio.*")

WORKING = $(shell pwd)
OBJDIR = .obj
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SRCS)))
DEPS = $(patsubst %.cpp,$(OBJDIR)/%.d,$(notdir $(SRCS)))
Dict = $(CODE)/ReadUrQMD.h
DictSrc = Dictionary.cxx

CXX = g++
OPTION = -O3 -std=c++11
DEBUG = -g -Wall
DEPFLAG = -MMD -MP
CC = $(CXX) ${OPTION} $(ROOTFLAG) $(DEBUG)

.PHONY: all clean
all: $(WORKING)/$(OBJDIR) run

run: $(OBJS) $(DictSrc)
	$(CC) $(ROOTLIBS) -o $@ $^

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CC) $(DEPFLAG) -c $< -o $@

$(DictSrc): $(Dict)
	rootcint -f $@ -c $^

$(WORKING)/$(OBJDIR):
	mkdir -p $(OBJDIR)

-include $(DEPS)

clean:
	rm -rf $(OBJDIR) Dict*
