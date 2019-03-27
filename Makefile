#compilers
CC = gcc
NVCC = nvcc

#flags 
CFLAGS = -Wall
<<<<<<< HEAD
NVCCFLAGS = -g -G --use_fast_math
=======
NVCCFLAGS = -g -G --use_fast_math -arch=sm_30
>>>>>>> c917ee52c29b4f0cf3f5a8867f0280edf316a78d

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
<<<<<<< HEAD
	$(RM) $(OBJECTS) $(CU_OBJECTS) $(TARGET) *.csv
=======
	$(RM) $(OBJECTS) $(CU_OBJECTS) $(TARGET) *.csv a.out
>>>>>>> c917ee52c29b4f0cf3f5a8867f0280edf316a78d
