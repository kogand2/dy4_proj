# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3
# CFLAGS = -c
LDFLAGS = -pthread
INCLUDE = -I../include/
SOURCES = threaded_stereo_rds.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp radio_fn.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = threaded_stereo_rds

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
