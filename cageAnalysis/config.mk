# This a common configuration file that includes definitions used by all
# the Makefiles.

# C++ compiler
CXX=g++

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -g -DNDEBUG -pipe
#CFLAGS=-O3 -Wall -ansi -pedantic -DNDEBUG

# Relative include and library paths for compilation of the examples
#E_INC=-I/Users/ismael.gomez/Documents/Software/Voro++/voro++-0.4.6/src
#E_LIB=-L/Users/ismael.gomez/Documents/Software/Voro++/voro++-0.4.6/src
#E_LIB= -lCGAL -lgmp
E_LIB= -lCGAL -lgmp -Lvoro/src -lvoro++

# Installation directory
PREFIX=/usr/local

# Install command
INSTALL=install

# Flags for install command for executable
IFLAGS_EXEC=-m 0755

# Flags for install command for non-executable files
IFLAGS=-m 0644
