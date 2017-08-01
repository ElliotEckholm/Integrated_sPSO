gcc -c SPSO.c -I/usr/local/include
gcc -c tools.c -I/usr/local/include
gcc -c mersenne.c -I/usr/local/include
gcc -c KISS.c -I/usr/local/include
gcc -c alea.c -I/usr/local/include
gcc -c wrapper.c -I/usr/local/include
gcc -c ptapso.c -I/usr/local/include
gcc -c maxphaseutils.c -I/usr/local/include
gcc -c ptapsotestfunc.c -I/usr/local/include
gcc -c test_ptapso.c -I/usr/local/include
gcc SPSO.o tools.o mersenne.o KISS.o alea.o wrapper.o ptapso.o maxphaseutils.o ptapsotestfunc.o test_ptapso.o -L/usr/local/lib -lgsl -lgslcblas -lm -o test_ptapso.exe
./test_ptapso.exe
