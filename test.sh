#!/bin/bash

rm -f *.csv

for i in 100 1000 5000 15360; do
	echo "Testing MxN: $i x $i, block_size: 16..."
    echo "/*--------------------------------------------*/"
    ./prog -n $i -m $i -p 1000 -a -s -g -t -w > a.out        
    for j in 32 64 128 256 1024; do
		echo "Testing MxN: $i x $i, block_size: $j..."
		echo "/*--------------------------------------------*/"
        ./prog -n $i -m $i -b $j -p 1000 -a -g -t -w > a.out
    done
done
