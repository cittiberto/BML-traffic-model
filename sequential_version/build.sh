#!/bin/bash
cd src
mkdir bld
cd bld
cmake -O3 ..
make
