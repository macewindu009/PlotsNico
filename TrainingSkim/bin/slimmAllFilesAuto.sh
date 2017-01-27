#specify directory name in which the skimmed data lies
gridname=/storage/jbod/nzaeh/2017-01-24/Gridoutput/

filenameTot=$( find $gridname -type d -name '*0000*' | grep 100to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0 && (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 2
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 200to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&&Boson_Pt > 70&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 3
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 400to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&&Boson_Pt > 90&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 4
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 600to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&&Boson_Pt > 110&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 5
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep M-50_T | grep amc)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 8
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep M-50_T | grep mad)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 1
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep M-150_ )
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 7
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 800to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& Boson_Pt > 130&&(recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 10 
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 1200to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& Boson_Pt > 150&&(recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 11
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep 2500to)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& Boson_Pt > 170&&(recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 12
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep M-10to | grep mad)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 6
	cp slimmed.root $filename
	rm slimmed.root
done

filenameTot=$( find $gridname -type d -name '*0000*' | grep M-10to | grep amc)
for filename in $filenameTot; do
	hadd $filename/outputtot.root -f $filename/output_*
	shrink $filename/outputtot.root "select > 0&& (recoilslimmedMETs_LongZ > -2000 && recoilslimmedMETs_LongZ < 2000)" MAPAnalyzer/t 9
	cp slimmed.root $filename
	rm slimmed.root
done

names=$(find $gridname -name "slimmed.root")
hadd -f $gridname/shrinkedDataTot.root $names

