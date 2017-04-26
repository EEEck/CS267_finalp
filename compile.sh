#for workstation
icpc --std=c++11 -fopenmp -o parallel_test_improved parallel_test_improved.cpp -O3 -larmadillo -g -parallel-source-info
#add for vector information
-qopt-report=1 -qopt-report-phase=vec -vec-report3
#for knl
icc  test.cpp  -o example -larmadillo -L/global/homes/m/mloipers/armadillo-7.800.2/build/usr/local/lib64 -I/global/homes/m/mloipers/armadillo-7.800.2/build/usr/local/include
