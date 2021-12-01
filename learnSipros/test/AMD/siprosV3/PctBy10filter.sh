#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate openmpi
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutBy10Pct & \
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutBy10Pct & \
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutBy10Pct
wait
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutBy10Pct & \
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutBy10Pct & \
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutBy10Pct
wait
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutBy10Pct & \
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutBy10Pct & \
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutBy10Pct