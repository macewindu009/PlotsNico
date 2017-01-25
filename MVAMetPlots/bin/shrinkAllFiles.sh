mergedDir=/storage/jbod/nzaeh/2017-01-17/Gridoutput/merged/
if [ ! -d "$mergedDir/shrinked/" ]; then
  # Control will enter here if $outputDir doesn't exist.
        mkdir $mergedDir/shrinked/
fi

for name in $( find $mergedDir -name '*.root'); do
	python shrink.py -i $name -o $mergedDir/shrinked/ -c mm mt &
done 
