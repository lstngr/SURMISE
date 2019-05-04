TARGET=surmise

CXX=g++ -std=c++11
CXXFLAGS=-g -Wall
LDFLAGS=
LDLIBS=

SRCDIR=src
INCDIR=include
OBJDIR=build
BINDIR=bin

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(INCDIR)/*.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

RM=rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

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
