#!/bin/bash

# Check if the docs directory exists, if not, create it
if [ ! -d "docs" ]; then
    mkdir docs
fi

# Compile the timing.c program
gcc -O3 -I./include src/matrix.c timing.c -o time_matvec -lm

# Run the program with different matrix sizes and save the output to timing.txt
{
    echo "Timing for m=n=1000"
    ./time_matvec 1000 1000 10

    echo "Timing for m=n=2000"
    ./time_matvec 2000 2000 10

    echo "Timing for m=n=3000"
    ./time_matvec 3000 3000 10

    echo "Timing for m=n=4000"
    ./time_matvec 4000 4000 10
} > docs/timing.txt
