
# complile openmp edition

VS Code C_C++ 项目快速配置模板   
- D:\微生物2021\2021-08\VS Code C_C++ 项目快速配置模板 - PRIN BLOG.mhtml
- https://printempw.github.io/vscode-cpp-project-quick-setup/

(By convention, gcc and other *nix compilers automatically add the lib prefix and the appropriate suffix.)
- https://stackoverflow.com/questions/12187078/why-cant-gcc-find-my-static-library

Using -no-pie should create a "regular" (position dependent) executable.
- https://stackoverflow.com/questions/46827433/g-compile-error-rodata-can-not-be-used-when-making-a-shared-object-recomp

```
sudo apt install openmpi-bin
```

make openmp use main.cpp  
make mpi use mpimain.cpp   
So I delete mpimain.cpp here

# test on server

solve undefined reference to `aligned_alloc@GLIBC_2.16' in conda
- https://stackoverflow.com/questions/68453991/gcc-linux-64-compiler-version-8-4-on-rhel-7-undefined-symbol-clock-gettimegli

```{bash}
cd /scratch/yixiong/learnSipros/
cp ../COVID19/MS2/MS2/Pan_050721_COVID1.ms2 MS2/
cp /scratch/yixiong/COVID19/maxquantOutPut20210716/monkeyCOVIDdb.fasta \
    /scratch/yixiong/COVID19/SiprosConfig.cfg \
    test
conda activate openmpi
cd /scratch/yixiong/Sipros-Ensemble/Scripts/
python sipros_prepare_protein_database.py \
  -i /scratch/yixiong/learnSipros/test/monkeyCOVIDdb.fasta \
  -o /scratch/yixiong/learnSipros/test/monkeyCOVIDdbDecoy.fasta \
  -c /scratch/yixiong/learnSipros/test/SiprosConfig.cfg

cd /scratch/yixiong/learnSipros/Sipros
mpic++ \
    -g \
    -std=gnu++11 \
    -o build/SiprosLearn \
    src/*.cpp \
    src/Scores/*.cpp \
    -Llib/ \
    -lmstoolkitlite \
    -lrt \
    -Wl,--sysroot=/ \
    -Isrc/MSToolkit/include/ \
    -DGCC \
    -D_FILE_OFFSET_BITS=64 \
    -fpermissive \
    -fopenmp \
    -O3 \
    -no-pie

cd /scratch/yixiong/learnSipros/test
mkdir output
nohup ../Sipros/build/SiprosLearn -o ./output \
  -w ../MS2 \
  -c SiprosConfig.cfg \
  > Sipros.log.txt 2>&1 &  
```

using openmp with cmake
  -https://stackoverflow.com/questions/54934933/undefined-reference-to-omp-get-wtime

cmake tutorial
  -https://www.jianshu.com/p/3078a4a195df
  -https://mp.weixin.qq.com/s/jgO7C7TwJ-yoBaxUQlXfyw
CMake 工程调用 Makefile 编译项目
  -https://blog.csdn.net/damenhanter/article/details/82499017