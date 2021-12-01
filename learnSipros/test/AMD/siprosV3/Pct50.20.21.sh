#!/bin/bash
for i in {20,21}; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel50 \
   -o MSlabel50Out
done
