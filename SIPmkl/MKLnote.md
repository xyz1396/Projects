
# documents

[MKL document](https://software.intel.com/content/www/us/en/develop/documentation/onemkl-linux-developer-guide/top/getting-started.html)  
[MKL document one api](https://docs.oneapi.io/versions/latest/onemkl/index.html)  
[conda](https://anaconda.org/anaconda/mkl)  
[install mkl on ubuntu](https://github.com/eddelbuettel/mkl4deb)  
[isotope abundance table](https://periodictable.com/Properties/A/IsotopeAbundances.html)

# install mkl by pypi

```bash
pip install mkl
pip install mkl-include
# or
pip install mkl-devel
# search funtion in .a or .so lib
cd /home/xyz/.local/lib
nm -A *.a | egrep '[TWBD] (dger_|dgemv_|dgetrf_|dgetri_)$'
nm -A *.so.* | egrep '[TWBD] (dger_|dgemv_|dgetrf_|dgetri_)$'
```

head file path: /home/xyz/.local/include/  
lib path: /home/xyz/.local/lib

[install mkl by conda](https://stackoverflow.com/questions/51109662/mkl-service-h-header-missing-from-anaconda-mkl-install)  
[install mkl by pip Intel one api document](https://software.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/pip.html)  
[search funtion in .a or .so lib](https://stackoverflow.com/questions/15935866/linking-intels-mkl-blas-lapack-to-gcc)  
[max or min of vector](https://stackoverflow.com/questions/9874802/how-can-i-get-the-maximum-or-minimum-value-in-a-vector)
[vector to pointer](https://stackoverflow.com/questions/6485496/how-to-get-stdvector-pointer-to-the-raw-data)

# demo
[MKL学习——基本操作C++实现](https://blog.csdn.net/zb1165048017/article/details/70216994)  
[Eigen库学习教程](https://blog.csdn.net/hongge_smile/article/details/107296658)  
  
[理解卷积](https://www.zhihu.com/question/22298352)  
[mkl convolution demo](https://gist.github.com/coder1234/4a1f661f1272dd6083e2)  
[Intel mkl convolution demo](https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/statistical-functions/convolution-and-correlation/convolution-and-correlation-usage-examples.html)