#!/bin/sh
../../akmc.py --reset -f
mpirun -n 8 ../../tools/emt_sp.py : \
       -n 8 ../../client/client : \
       -n 1 ../../akmc.py 
