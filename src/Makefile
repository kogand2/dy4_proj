# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3
# CFLAGS = -c
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES = fmRDSDesign.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp block_conv_fn.cpp mono_path.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = fmRDSDesign

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
