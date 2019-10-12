mpirun --oversubscribe -np 8 ./build/Cartint 2 2 1.0 1.0 > w1_d2_L2.txt
mpirun --oversubscribe -np 8 ./build/Cartint 2 6 1.0 1.0 > w1_d2_L6.txt
mpirun --oversubscribe -np 8 ./build/Cartint 2 12 1.0 1.0 > w1_d2_L12.txt
mpirun --oversubscribe -np 8 ./build/Cartint 2 20 1.0 1.0 > w1_d2_L20.txt
mpirun --oversubscribe -np 8 ./build/Cartint 2 30 1.0 1.0 > w1_d2_L30.txt
mpirun --oversubscribe -np 8 ./build/Cartint 2 42 1.0 1.0 > w1_d2_L42.txt
