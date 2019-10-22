#!/bin/bash
gcc -o main main.c -fopenmp -O3 -I../mylib/include/ ../mylib/obj/*.o -lm rust_code2/target/release/librl_pendulum.so
