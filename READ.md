# Integrated_sPSO


To run the code: 


gcc -c ptapso.c -I/usr/local/include
gcc -c maxphaseutils.c -I/usr/local/include
gcc -c ptapsotestfunc.c -I/usr/local/include
gcc -c test_ptapso.c -I/usr/local/include
gcc ptapso.o maxphaseutils.o ptapsotestfunc.o test_ptapso.o -L/usr/local/lib -lgsl -lgslcblas -lm -o test_ptapso.exe
./test_ptapso.exe
