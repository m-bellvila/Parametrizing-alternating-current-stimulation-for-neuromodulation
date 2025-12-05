This code was used in the data analysis for:

Monica Bell Vila, Margaret Koletar, Adrienne Dorr, Andrea Trevisiol, Maged Goubran, John G. Sled, Paolo Bazzigaluppi, Bojana Stefanovic; 
Parametrizing alternating current stimulation for neuromodulation. Imaging Neuroscience 2025; doi: https://doi.org/10.1162/IMAG.a.1066


MUA analysis:
First the sACS artefact was removed from MUA data using MUAartrem.m 
Data was concatenated for spike detection using concatdataforspikesorting.py 
After discarding overly channels or stimulation periods with incomplete artefact removal, data was submitted to SpyKING CIRCUS 



Following spike detection & neuron clustering using SpyKING CIRCUS:
phaseextract.py 
nulldistribution.py
signeurons.py 
