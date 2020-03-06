#!/bin/bash
rm obj/*.o
echo "Compiling the source files in src/ ..." &&
#gcc -c -g -Wextra -Wshadow -pedantic -Wno-strict-overflow -Wall -Iinclude/ src/*.c &&

# Compiling option extremely pedantic:
#gcc -c -fopenmp -Wextra -Wshadow -pedantic -Wall -O3 -Iinclude/ src/*.c &&

# Compiling option fine for my usage, where I do not care about
# shadowed variables, unused values or optimization warning
gcc -c -fopenmp -Wextra -Wno-unused-value -pedantic -Wall -Wno-strict-overflow-O3 -Iinclude/ src/*.c &&
echo "Moving the objects into obj/ ..." &&
mv *.o obj/ &&
echo "Done!"
