#!/bin/bash
rm obj/*
echo "Compiling the source files in src/ ..." &&
#gcc -c -g -Wextra -Wshadow -pedantic -Wall -Iinclude/ src/*.c &&
gcc -c -fopenmp -Wextra -Wshadow -pedantic -Wall -O3 -Iinclude/ src/*.c &&
echo "Moving the objects into obj/ ..." &&
mv *.o obj/ &&
echo "Done!"
