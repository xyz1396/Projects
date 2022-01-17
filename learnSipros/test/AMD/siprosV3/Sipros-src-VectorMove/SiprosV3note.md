[mpi omp混合编程](https://www.lainme.com/doku.php/blog/2014/12/mpi%E5%92%8Copenmp%E6%B7%B7%E5%90%88%E7%BC%96%E7%A8%8B%E7%AC%94%E8%AE%B0)  
[slurm tutorial](https://bicmr.pku.edu.cn/~wenzw/pages/slurm.html)  


# local compile

[mpich portable mpi](https://www.cnblogs.com/xingkongyihao/p/9733260.html)

```bash
sudo apt install mpich
cd build
cmake ..
make -j6
```

# on server

```bash
conda install -c bioconda biopython=1.68

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3

python scripts/reverseseq.py \
  -i /scratch/yixiong/learnSipros/test/AMD/fungiNoDup.faa \
  -o fungiDecoy.faa \

python scripts/sip.py \
    -c SiprosConfig.N15_SIP.cfg \
    -w configFiles

nohup ./SiprosV3omp \
  -g configFiles \
  -w ../MSlabel0 \
  -o MSlabel0Out \
  > MSlabel0Out.log.txt 2>&1 & 

#!/bin/bash
#SBATCH --partition=omicsbio
#SBATCH --ntasks=60
#SBATCH --time=999:00:00
#SBATCH --job-name=MSlabel98
#SBATCH --mem=150G
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD/siprosV3
#SBATCH --output=MSlabel98Out_stdout.txt
#SBATCH --error=MSlabel98Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3
./Sipros \
  -g configFiles \
  -w ../MSlabel98 \
  -o MSlabel98Out

vim MSlabel98Out.sb

./SiprosV3omp -c configurefilename -w workingdirectory

#!/bin/bash
export OMP_NUM_THREADS=40
for i in `seq 0 1 100`; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel0 \
   -o MSlabel0Out
done

vim Pct0.sh
nohup ./Pct0.sh > MSlabel0Out.log.txt 2>&1 &

#!/bin/bash
export OMP_NUM_THREADS=40
for i in `seq 0 1 100`; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel50 \
   -o MSlabel50Out
done

vim Pct50.sh
nohup ./Pct50.sh > MSlabel50Out.log.txt 2>&1 & 

#!/bin/bash
for i in {20,21}; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel50 \
   -o MSlabel50Out
done

vim Pct50.20.21.sh
nohup ./Pct50.20.21.sh > MSlabel50Out20.21.log.txt 2>&1 &

#!/bin/bash
export OMP_NUM_THREADS=40
for i in `seq 0 1 100`; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel98 \
   -o MSlabel98Out
done

vim Pct98.sh
nohup ./Pct98.sh > MSlabel98Out.log.txt 2>&1 & 

#!/bin/bash
for i in {20,21}; do
   ./SiprosV3omp -c configFiles/SiproConfig.N15_${i}Pct.cfg \
   -w ../MSlabel98 \
   -o MSlabel98Out
done

vim Pct98.20.21.sh
nohup ./Pct98.20.21.sh > MSlabel98Out20.21.log.txt 2>&1 &

#!/bin/bash
#SBATCH --time=999:00:00
#SBATCH --job-name=MSlabel50
#SBATCH --partition=omicsbio
#SBATCH --nodes=2 
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=10
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD/siprosV3
#SBATCH --output=MSlabel50Out_stdout.txt
#SBATCH --error=MSlabel50Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=10
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3
mpirun -np 8 ./SiprosV3mpi \
  -g configFiles \
  -w ../MSlabel50 \
  -o MSlabel50Out

vim MSlabel50Out.sb  

#!/bin/bash
#SBATCH --time=999:00:00
#SBATCH --job-name=MSlabel50
#SBATCH --partition=omicsbio
#SBATCH --nodes=2 
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=10
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD/siprosV3
#SBATCH --output=MSlabel50Out_stdout.txt
#SBATCH --error=MSlabel50Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

module load OpenMPI/4.0.3-GCC-9.3.0
export OMP_NUM_THREADS=10
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3
mpirun -np 8 ./SiprosV3mpi \
  -g configFiles \
  -w ../MSlabel50 \
  -o MSlabel50OutMPI

vim MSlabel50Out.sb  

#!/bin/bash
#SBATCH --time=999:00:00
#SBATCH --job-name=MSlabel50
#SBATCH --partition=omicsbio
#SBATCH --nodes=2 
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=40
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD/siprosV3
#SBATCH --output=MSlabel50Out_stdout.txt
#SBATCH --error=MSlabel50Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=40
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3
mpirun ./SiprosV3mpi \
  -g configFiles \
  -w ../MSlabel50 \
  -o MSlabel50Out
```

# Performance test

all cores are running   
less cores will get higher single cpu utilization

```bash
ssh c659
export OMP_NUM_THREADS=5
./SiprosV3ompTest -c configFiles/SiproConfig.N15_0Pct.cfg    -w ../MSlabel98    -o /dev/null -s
export OMP_NUM_THREADS=40
./SiprosV3ompTest -c configFiles/SiproConfig.N15_0Pct.cfg    -w ../MSlabel98    -o /dev/null -s
```

# Increase core efficiency

change PEPTIDE_ARRAY_SIZE from 10000 to 200000

```bash
cd /scratch/yixiong/learnSipros/test/AMD/siprosV3/Sipros-src
# 4398% cpu
bin/SiprosV3omp1000 -c ../configFiles/SiproConfig.N15_0Pct.cfg -w ../../MSlabel98 -o /dev/null -s
# 1764% cpu
bin/SiprosV3omp10000 -c ../configFiles/SiproConfig.N15_0Pct.cfg -w ../../MSlabel98 -o /dev/null -s
# 1000-1500% cpu
bin/SiprosV3omp200000 -c ../configFiles/SiproConfig.N15_0Pct.cfg -w ../../MSlabel98 -o /dev/null -s
```

## test real run time

```bash
cd /scratch/yixiong/learnSipros/test/AMD/siprosV3/Sipros-src/timeCompare

#!/bin/bash
echo "1k"
time ./SiprosV3omp1000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 1k -s
echo "10k"
time ./SiprosV3omp10000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 10k -s
echo "200k"
time ./SiprosV3omp200000 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -o 200k -s

vim compareTime.sh
nohup ./compareTime.sh > compareTime.log.txt 2>&1 &
```

1k
real    2m38.849s
user    113m56.614s
sys     0m1.352s

10k
real    2m9.283s
user    36m33.607s
sys     0m1.939s

200k
real    2m4.954s
user    26m27.802s
sys     0m7.460s

## gperftools

```bash
# apt install 
sudo apt install google-perftools libgoogle-perftools-dev

module gperftools/2.7.90-GCCcore-8.3.0

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3/Sipros-src/timeCompare
CPUPROFILE=./Gprofile
./SiprosV3ompGper -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s
# use browser
pprof SiprosV3ompGper  test_capture.prof --web
pprof SiprosV3ompGper test_capture.prof --pdf > prof.pdf
pprof SiprosV3ompGper test_capture.prof --text > prof.txt
```
## optimise vector

[vector optimise](https://stackoverflow.com/questions/49615076/is-the-poor-performance-of-stdvector-due-to-not-calling-realloc-a-logarithmic)

```bash
# real    2m56.632s
# user    53m3.296s
# sys     0m7.480s
time ./SiprosV3ompOldSum -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    2m29.344s
# user    37m48.363s
# sys     0m7.450s
time ./SiprosV3ompOldSumReserve -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    2m36.476s
# user    51m43.590s
# sys     0m5.889s
time ./SiprosV3ompNewSum -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    1m58.424s
# user    32m21.408s
# sys     0m5.609s
time ./SiprosV3ompNewSumReserve -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    1m57.984s
# user    31m56.847s
# sys     0m5.924s
time ./SiprosV3ompNewSumStatic -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    2m1.196s
# user    32m47.772s
# sys     0m5.934s
time ./SiprosV3ompGperO3 -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    2m26.035s
# user    46m13.126s
# sys     0m6.632s
time ./SiprosV3ompGperNoReserve -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    1m58.517s
# user    31m30.472s
# sys     0m5.834s
time ./SiprosV3ompReserve -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# real    2m20.190s
# user    44m6.703s
# sys     0m5.787s
time ./SiprosV3ompNoReserve -c SiproConfig.N15_0Pct.cfg -f AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 -s

# running time： 126s
# vector
./make run

# running time： 131s
# vector with move
./make run
```

# Filter

```bash
#!/bin/bash
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0Out

python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0Out

python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0Out

vim Pct0filter.sh
nohup ./Pct0filter.sh > Pct0filter.log.txt 2>&1 & 
```

# Filter step by 10

```bash
cd /scratch/yixiong/learnSipros/test/AMD/siprosV3/
for i in `seq 0 10 100`; do
   cp MSlabel0Out/*_${i}Pct* MSlabel0OutBy10Pct
   cp MSlabel50Out/*_${i}Pct* MSlabel50OutBy10Pct
   cp MSlabel98Out/*_${i}Pct* MSlabel98OutBy10Pct
done

for i in `seq 5 10 100`; do
   cp MSlabel0Out/*_${i}Pct* MSlabel0OutFrom5Pct
   cp MSlabel50Out/*_${i}Pct* MSlabel50OutFrom5Pct
   cp MSlabel98Out/*_${i}Pct* MSlabel98OutFrom5Pct
done

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

vim PctBy10filter.sh
nohup ./PctBy10filter.sh > PctBy10filter.log.txt 2>&1 &

#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate openmpi
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutFrom5Pct & \
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutFrom5Pct & \
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutFrom5Pct
wait
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutFrom5Pct & \
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutFrom5Pct & \
python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutFrom5Pct
wait
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel0OutFrom5Pct & \
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50OutFrom5Pct & \
python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel98OutFrom5Pct

vim PctFrom5filter.sh
nohup ./PctFrom5filter.sh > PctFrom5filter.log.txt 2>&1 &
```