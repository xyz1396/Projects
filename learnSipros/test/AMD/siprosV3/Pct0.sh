#!/bin/bash
export OMP_NUM_THREADS=40
for i in `seq 0 1 100`; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel0 \
   -o MSlabel0Out
done
