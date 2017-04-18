mergedDir=/nfs/dust/cms/user/nzaeh/2016-12-13_11-23_analysis/merged/
outputDir=$mergedDir/CSVFiles
if [ ! -d "$outputDir" ]; then
  # Control will enter here if $outputDir doesn't exist.
	mkdir $outputDir
fi
for filename in $( find $mergedDir -name '*.root'); do
	makeCSVFromArtus $filename $outputDir mm mt
done
