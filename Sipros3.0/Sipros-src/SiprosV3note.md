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

module load CMake/3.15.3-GCCcore-8.3.0\
  OpenMPI/3.1.4-GCC-8.3.0                        
export OMP_NUM_THREADS=10
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

cd /scratch/yixiong/learnSipros/test/AMD/siprosV3/
mpirun -np 8 ./Sipros-src/bin/SiprosV3mpi \
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

#!/bin/bash
python scripts/sipros_peptides_filtering.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50Out

python scripts/sipros_peptides_assembling.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50Out

python scripts/ClusterSip.py \
  -c SiprosConfig.N15_SIP.cfg \
  -w MSlabel50Out

vim Pct50filter.sh
nohup ./Pct50filter.sh > Pct50filter.log.txt 2>&1 &

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

vim Pct98filter.sh
nohup ./Pct98filter.sh > Pct98filter.log.txt 2>&1 &
```