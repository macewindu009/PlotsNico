mergedDir=/storage/jbod/nzaeh/2017-01-24/
for name in $( find $mergedDir -name '*24_0*'); do
	rm -r $name 
done 

#for name in $( find /storage/jbod/nzaeh/2017-01-24/ -name '*.root'); do
       #printf $name 
#done 
