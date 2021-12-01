#!/bin/bash
echo "1k"
time ./SiprosV3omp1000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 1k -s
echo "10k"
time ./SiprosV3omp10000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 10k -s
echo "200k"
time ./SiprosV3omp200000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 200k -s
