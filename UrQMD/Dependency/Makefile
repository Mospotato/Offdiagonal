ROOTFLAG = `root-config --cflags`
ROOTLIBS = `root-config --libs`

LDFLAGS = -fPIC
OBJDIR = .obj
CXX = g++ -O3 -std=c++11
DEBUG = -Wall -pedantic
DEPFLAG = -MMD -MP

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(notdir $(SRCS)))
DEPS = $(patsubst %.cpp,$(OBJDIR)/%.d,$(notdir $(SRCS)))
CC = $(CXX) $(DEBUG) $(DEPFLAG) $(INCLUDE)
LIBRARY = libDependency.so
.PHONY: all clean
all: $(LIBRARY)
$(LIBRARY): $(OBJS) 
	$(CC) $(LDFLAGS) $(ROOTLIBS) -shared -o $@ $^ 

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CC) $(ROOTFLAG) $(LDFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(CONFIG)/%.cpp | $(OBJDIR)
	$(CC) $(ROOTFLAG) $(LDFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

-include $(DEPS)

clean:
	rm -rf $(OBJDIR)
