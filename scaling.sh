#!bin/bash

export OMP_NUM_THREADS=1
./parallel_test > slow.out
./parallel_test_improved >fast.out

export OMP_NUM_THREADS=2
./parallel_test >> slow.out
./parallel_test_improved >>fast.out

export OMP_NUM_THREADS=4
./parallel_test >> slow.out
./parallel_test_improved >>fast.out

export OMP_NUM_THREADS=8
./parallel_test >> slow.out
./parallel_test_improved >>fast.out

export OMP_NUM_THREADS=16
./parallel_test >> slow.out
./parallel_test_improved >>fast.out

export OMP_NUM_THREADS=32
./parallel_test >> slow.out
./parallel_test_improved >>fast.out
