libTHTensor.0.so:
	gcc -O3 -msse4.2 -mfma -mavx -mavx2 -march=native -mtune=native -fPIC -fopenmp -c *.c vector/AVX.c vector/AVX2.c
	gcc -fopenmp -shared -o libTHTensor.0.so *.o -lopenblas -lrt
