libTHNeural.0.dylib:
	gcc -O3 -msse4.2 -mfma -mavx -mavx2 -march=native -mtune=native -c -I../TH -I. init.c
	gcc -shared -o libTHNeural.0.dylib *.o -L/usr/local/lib -lTHTensor
