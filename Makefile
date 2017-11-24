# ------------------------------------------------
# Generic Makefile
# ------------------------------------------------

# project name (generate executable with this name)
TARGET   = torsionEvol

CC	   = g++
# compiling flags here
CFLAGS   = -I. -O3 -g
CXXFLAGS = -std=c++11 -Wall -Wpedantic -Wextra

LINKER   = g++
# linking flags here
LDFLAGS   = -I. -lm -ldl -lboost_system -lboost_filesystem 

# change these to proper directories where each file should be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm	   = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(LINKER) $(OBJECTS) $(LDFLAGS) -o $@
	@echo "\033[00;32mLinking completed.\033[00m"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(OBJDIR)
	$(CC) -c $(CFLAGS) $(CXXFLAGS) $< -o $@
	@echo "\033[00;32mCompiled "$<".\033[00m"

.PHONY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "\033[00;32mCleanup completed.\033[00m"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "\033[00;32mExecutable removed!\033[00m"