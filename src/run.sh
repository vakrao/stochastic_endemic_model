gcc -c *.c 
echo '#### Compiling model! ####' 
gcc -o model *.o -lgsl -lgslcblas -lm 
echo '#### Running model! ####' 
./model 1 00 0 one_strain.csv ../data/
