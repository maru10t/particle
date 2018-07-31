#!/bin/sh
g++ -std=c++11 --shared -fPIC -O2 -shared *.cc -o libmd.so -fopenmp
