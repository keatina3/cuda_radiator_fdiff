#!/bin/bash

rm -f *.csv

<<<<<<< HEAD
for i in 100 1000 5000 15360; do       
=======
for i in 100 1000 5000 15360; do
	#echo "Testing MxN: $i x $i, block_size: 16..."
	#echo "/*--------------------------------------------*/"
    #./prog -n $i -m $i -p 1000 -a -s -g -t -w > a.out        
>>>>>>> c917ee52c29b4f0cf3f5a8867f0280edf316a78d
    for j in 16 32 64 128 256 1024; do
		echo "Testing MxN: $i x $i, block_size: $j..."
		echo "/*--------------------------------------------*/"
        ./prog -n $i -m $i -b $j -p 1000 -a -d -t -w > a.out
    done
done
