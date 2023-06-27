gcc -c *.c 
echo '#### Compiling model! ####' 
gcc -o model *.o -lgsl -lgslcblas -lm 
#echo '#### Running model! ####' 
./model 1 30 "category" one_strain.csv ../ "onerun"
