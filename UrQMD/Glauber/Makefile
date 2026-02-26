### ROOT setting
ROOTLIBS = `root-config --libs`
ROOTFLAG = `root-config --cflags`
### Configuration
CODE = Code
OBJDIR = $(CODE)/obj
SRCS = $(wildcard $(CODE)/*.cxx)
OBJS = $(patsubst %.cxx,$(OBJDIR)/%.o,$(notdir $(SRCS)))
DEPS = $(patsubst %.cxx,$(OBJDIR)/%.d,$(notdir $(SRCS)))
### FLAGS
CXX = g++
OPTION = -O3
DEBUG = -g -Wall -Wextra -pedantic
DEPFLAG = -MMD -MP
INCLUDE = -I$(CODE)
CC = $(CXX) ${OPTION} $(ROOTFLAG) $(INCLUDE)
### targets
all: run

run: $(OBJS)
	$(CC) $(ROOTLIBS) -o $@ $^

$(OBJDIR)/%.o: $(CODE)/%.cxx | $(OBJDIR)
	$(CC) $(DEPFLAG) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

-include $(DEPS)

### clean up
clean:
	rm -rf $(OBJDIR)/*
	ls
