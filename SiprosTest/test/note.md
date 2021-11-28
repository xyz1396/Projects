

[Sipros V3 code](https://code.google.com/archive/p/sipros/source/default/source)  
[Sipros V3 tutorial](https://sipros.omicsbio.org/instructions/sip_instructions/)

# SIP test

```bash
cd /scratch/yixiong/Sipros-Ensemble/Scripts/
python sipros_prepare_protein_database.py \
  -i /scratch/yixiong/learnSipros/test/Ecoli/Ecoli.fasta \
  -o /scratch/yixiong/learnSipros/test/Ecoli/EcoliDecoy.fasta \
  -c /scratch/yixiong/learnSipros/test/Ecoli/SiprosConfigC13SIP.cfg
# this will generate 100 cfg files at different isotopic percent
python sip.py \
    -c /scratch/yixiong/learnSipros/test/Ecoli/SiprosConfigC13SIP.cfg \
    -w /scratch/yixiong/learnSipros/test/Ecoli/
cd /scratch/yixiong/learnSipros/test
./Sipros -o ./output \
  -f ../SIP/Ecoli_OU.ms2 \
  -c SiprosConfigC13SIP.cfg
./Sipros -o ./output50 \
  -f ../SIP/Ecoli_OU.ms2 \
  -c SiproConfig.C13_50Pct.cfg
```

# SIP Experiment

AMD_0Percent15N_Velos_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210 [0% N15 velos]
AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110 [98% N15 velos]

```bash

# make database
cd /scratch/yixiong/learnSipros/test/AMD
cp /scratch/yixiong/learnSipros/SIP/StandardEnrichment/all_c75uniq_v4_Core_Fungal_FR.fasta .
seqkit grep all_c75uniq_v4_Core_Fungal_FR.fasta -r -p ^Rev_ -n -v -o fungi.faa
# 970 duplicated records removed
seqkit rmdup fungi.faa -s -o fungiNoDup.faa -D fungiDup.txt

cd /scratch/yixiong/Sipros-Ensemble/Scripts/
python sipros_prepare_protein_database.py \
  -i /scratch/yixiong/learnSipros/test/AMD/fungiNoDup.faa \
  -o /scratch/yixiong/learnSipros/test/AMD/fungiDecoy.faa \
  -c /scratch/yixiong/learnSipros/test/AMD/SiprosConfigN15SIP.cfg
python sip.py \
    -c /scratch/yixiong/learnSipros/test/AMD/SiprosConfigN15SIP.cfg \
    -w /scratch/yixiong/learnSipros/test/AMD/configFiles

cd /scratch/yixiong/learnSipros/test/AMD
cp /scratch/yixiong/learnSipros/SIP/StandardEnrichment2/velos_Orbi/AMD_0Percent15N_Velos_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210/*.FT2 MSlabel0
cp /scratch/yixiong/learnSipros/SIP/StandardEnrichment2/velos_Orbi/AMD_StandardEnrichment_SampleBeta_50Percent15N_Velos_OrbiMS2_Run2_020110/*.FT2 MSlabel50
cp /scratch/yixiong/learnSipros/SIP/StandardEnrichment2/velos_Orbi/AMD_StandardEnrichment_SampleAlpha_98Percent15N_Velos_OrbiMS2_Run2_013110/*.FT2 MSlabel98

nohup ../Sipros \
  -g configFiles \
  -w MSlabel0 \
  -o MSlabel0Out \
  > MSlabel0Out.log.txt 2>&1 & 
nohup ../Sipros \
  -g configFiles \
  -w MSlabel50 \
  -o MSlabel50Out \
  > MSlabel50Out.log.txt 2>&1 &
nohup ../Sipros \
  -g configFiles \
  -w MSlabel98 \
  -o MSlabel98Out \
  > MSlabel98Out.log.txt 2>&1 &

#!/bin/bash
#SBATCH --partition=omicsbio
#SBATCH --ntasks=40
#SBATCH --time=999:00:00
#SBATCH --job-name=MSlabel98
#SBATCH --mem=150G
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD
#SBATCH --output=MSlabel98Out_stdout.txt
#SBATCH --error=MSlabel98Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

cd /scratch/yixiong/learnSipros/test/AMD
../Sipros \
  -g configFiles \
  -w MSlabel98 \
  -o MSlabel98Out

vim MSlabel98Out.sb

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel0Out \
  -o MSlabel0OutFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel0OutFiltered.log.txt 2>&1 &

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel50Out \
  -o MSlabel50OutFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel50OutFiltered.log.txt 2>&1 &

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel98Out \
  -o MSlabel98OutFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel98OutFiltered.log.txt 2>&1 &
```

# step by 10 percent

```bash
for i in `seq 0 10 100`; do
   cp MSlabel0Out/*_${i}Pct* MSlabel0OutBy10Pct
done

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel0OutBy10Pct \
  -o MSlabel0OutBy10PctFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel0OutBy10PctFiltered.log.txt 2>&1 &

for i in `seq 0 10 100`; do
   cp MSlabel50Out/*_${i}Pct* MSlabel50OutBy10Pct
done

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel50OutBy10Pct \
  -o MSlabel50OutBy10PctFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel50OutBy10PctFiltered.log.txt 2>&1 &

for i in `seq 0 10 100`; do
   cp MSlabel98Out/*_${i}Pct* MSlabel98OutBy10Pct
done

nohup /scratch/yixiong/Sipros-Ensemble/Scripts/runSiprosFiltering.sh \
  -in MSlabel98OutBy10Pct \
  -o MSlabel98OutBy10PctFiltered \
  -c SiprosConfigN15SIP.cfg \
  > MSlabel98OutBy10PctFiltered.log.txt 2>&1 &
```