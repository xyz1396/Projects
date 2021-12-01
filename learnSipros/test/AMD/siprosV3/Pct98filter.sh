#!/bin/bash
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98Out

python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98Out

python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98Out
