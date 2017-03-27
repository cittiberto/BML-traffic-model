#!/bin/bash
cd src
mkdir bld
cd bld
cmake -fopenmp -O3 ..
make
