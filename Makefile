TARGET=surmise
export DOXY_VERSION=$(shell git describe --abbrev=8 --dirty --always --tags)
export OMPI_CXX=g++
#export OMPI_CXX=clang++

CXX=mpicxx
CXXFLAGS=-std=c++11 -g -Wall -pedantic -O3# -D_GLIBCXX_DEBUG #-fno-stack-protector
LDFLAGS=
LDLIBS=

SRCDIR=src
INCDIR=include
OBJDIR=build
BINDIR=bin
CXXFLAGS+= -I./$(INCDIR)

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(INCDIR)/*.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

RM=rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: tags
tags:
	$(RM) tags
	ctags $(SOURCES) $(INCLUDES)

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
	@(echo "Cleaned up objects.")

.PHONY: cleanall
cleanall: clean
	$(RM) $(BINDIR)/$(TARGET)
	@(echo "Cleaned up binaries.")

.PHONY: docs
docs:
	doxygen
