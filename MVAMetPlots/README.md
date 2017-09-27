# Creating Plots

Plots can be created with makePlotsMVAMet.py by using the optional arguments

It can for example be run with
 python makePlotsMVAMet.py -i /storage/jbod/nzaeh/2017-02-02/Gridoutput/WithEverything/dataMVAMetWithEverything.csv -o FinalTrainingMMWithPuppiQuantile --export -m Quantile -s -c '(1==plotData[:,dictPlot["select"]])' --fixedRange -d '$H \to \tau\tau \to \mu\tau_h$'

Fixed ranges and positioning of textboxes have to be adjusted in the code!

HiggsPlots.py does the same with replacing Z by H
