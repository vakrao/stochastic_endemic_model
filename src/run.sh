gcc -c *.c 
echo '#### Compiling model! ####' 
gcc -o model *.o -lgsl -lgslcblas -lm 
#echo '#### Running model! ####' 
#./model 1 30 2 one_strain.csv ../data/
