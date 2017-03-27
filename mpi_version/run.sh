#!/bin/bash
foo=$(nproc)
mpirun -n $foo ./out_MPI
