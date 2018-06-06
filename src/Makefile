# generally useful things

SHELL = /bin/sh

# ---------------------------------------------------------------------
# compiler/linker settings
# specify flags and libraries needed for your compiler

CC =		g++
CCFLAGS =	-Ofast
SHFLAGS =
DEPFLAGS =	-M

LINK =		g++
LINKFLAGS =	-Ofast
LIB =

ARCHIVE =	ar
ARFLAGS =	-rc
SHLIBFLAGS =	-shared
# ---------------------------------------------------------------------
# platform specific stuff

PLATFORM = $(shell uname)
ifeq ($(PLATFORM), Linux)
	EXTRA_COMPILE_FLAGS = -std=c++11 -DHAVE_INLINE
	LIBRARIES = -L/usr/lib -lgsl -lgslcblas -lm
else
	EXTRA_COMPILE_FLAGS = -DHAVE_INLINE -I/usr/local/include
	LIBRARIES = -L/usr/local/lib -lgsl -lgslcblas -lm
endif

CCFLAGS += $(EXTRA_COMPILE_FLAGS)


# ---------------------------------------------------------------------
# location and executable name settings and other such things...

BIN_DIR = ../bin
EXE = $(BIN_DIR)/vsim

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)

OBJ = $(SRC:.cpp=.o)

DOC_DIR = ../doc
DOC = $(DOC_DIR)/html

# ---------------------------------------------------------------------
# Link target

default: $(EXE)

clean:
	rm -rf *.o
	rm -rf *.d
	rm -rf $(EXE)

clean-doc:
	rm -rf $(DOC)

clean-all: clean clean-doc

doc:
	doxygen $(DOC_DIR)/Doxyfile

$(EXE):	$(OBJ) $(EXTRA_LINK_DEPENDS)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE) $(LIBRARIES)

# Compilation rules

%.o:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(SHFLAGS) $(EXTRA_INC) -c $< -o $@

%.d:%.cpp $(EXTRA_CPP_DEPENDS)
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)
