gcc -Wall -I/usr/local/include -c  -malign-double -fno-math-errno  gsl_test.c  
gcc -L/usr/local/lib gsl_test.o -lgsl -lgslcblas -lm -o gsl_test
