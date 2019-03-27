#compilers
CC = gcc
NVCC = nvcc

#flags 
CFLAGS = -Wall
NVCCFLAGS = -g -G --use_fast_math

#files
OBJECTS = main.o findiff.o radiator.o utils.o
CU_OBJECTS = findiff_gpu.o
CU_SOURCES = findiff_gpu.cu

TARGET = prog

all: $(OBJECTS) cu_objs
	$(NVCC) $(OBJECTS) $(CU_OBJECTS) -o $(TARGET)

cu_objs: $(CU_SOURCES)
	$(NVCC) $(CU_SOURCES) -c $(NVCCFLAGS)

.PHONY: clean

test: all
	./test.sh

clean:
	$(RM) $(OBJECTS) $(CU_OBJECTS) $(TARGET) *.csv a.out
