ROOTFLAG = `root-config --cflags`
ROOTLIBS = `root-config --libs`

CODE = Code
DEPDIR = ../Dependency
VPATH = $(CODE):$(DEPDIR)
SRCS = $(shell find $(CODE) $(DEPDIR) -name "*.cpp")

OBJDIR = $(shell pwd)/.obj
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SRCS)))
DEPS = $(patsubst %.cpp,$(OBJDIR)/%.d,$(notdir $(SRCS)))

CXX = g++
OPTION = -O3 -std=c++11
DEBUG = -g -Wall -Wextra
DEPFLAG = -MMD -MP
CC = $(CXX) ${OPTION} $(ROOTFLAG) $(DEBUG)

.PHONY: all clean
all: run

run: $(OBJS)
	$(CC) $(ROOTLIBS) -o $@ $^

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CC) $(DEPFLAG) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

-include $(DEPS)

clean:
	rm -rf $(OBJDIR)
	ls
