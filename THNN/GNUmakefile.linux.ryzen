libTHNeural.0.dylib:
	gcc -O3 -msse4.2 -mfma -mavx -mavx2 -march=native -mtune=native -fPIC -fopenmp -c -I../TH -I. init.c
	gcc -fopenmp -shared -o libTHNeural.0.so *.o -L/usr/local/lib -lTHTensor
