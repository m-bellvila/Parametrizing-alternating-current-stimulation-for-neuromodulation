#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#This script is to be run for data following processing of concatenated MUA
#through SpyKING Circus 
#It serves to extract the phases of each neuron cluster based on the spike 
#times relative to the stimulation waveform
#Only neurons with firing rates of at least 0.25Hz are stored for subsequent analysis


#%%
import numpy as np
import scipy.io
import h5py
import mat73
from scipy.signal import hilbert
import pickle
import glob



#%%
#Insert location of 1-min sACS experiment folders
folders = '/1minsACS/202*/'


Fs = 24414.0625 #Sampling frequency of MUA data

def phasescalc(stimwf,spiks):
    spikarr = np.copy(spiks)
    while len(np.where(spikarr >= len(stimwf))[0]) > 0:
        #If spike time is outside the duration of the stimulation waveform array, 
        #restart at beginning of stimulation waveform
        listtofix = np.where(spikarr >= len(stimwf))[0]
        print(len(stimwf))
        for bb in listtofix:
            print('Before: '+ str(spikarr[bb]))
            spikarr[bb] = spikarr[bb]- len(stimwf)
            print('After: ' +str(spikarr[bb]))
    #Hilbert transfrom to extract phase of stimulation
    dhilbert = hilbert(stimwf) #form the analytical signal
    stimphase = np.angle(dhilbert) #inst phase
    

    phasespike= np.zeros((len(spikarr))) #what is the stim phase at spike time
    ampspike = np.zeros((len(spikarr))) #what is the amplitude of stim at spike time

    j = []
    for j in range(len(spikarr)):
        spike = int(spikarr[j])
        phasespike[j]= stimphase[spike]  
        ampspike[j]= stimwf[spike]
    return [phasespike, ampspike]


for folder in folders:
    print(folder)
    segmentlengths = np.load(folder + 'segmentlengths.npy', allow_pickle = True)
    
    
    for i in range(len(segmentlengths[0])-1):
        segmentlengths[0][i+1] = segmentlengths[0][i+1] + segmentlengths[0][i]

    stimtp = segmentlengths[1]    
    amp = segmentlengths[4]
    freq = segmentlengths[3]
    stimwaveform = segmentlengths[2]
    
    
    clusterfile = folder + 'MUAfiltdata/MUAfiltdata.result-merged.hdf5'
    c = h5py.File(clusterfile, 'r')['spiketimes'] 
    clusters = list(c.keys())
    

    allphases = []
    storespikes = {}
    
    i = []
    tpcounter = np.zeros(len(stimtp))
    for i in range(len(clusters)):
        print('Cluster: '+str(i))
        spiketimes = np.array(sorted(c[clusters[i]]))
        clusterlabel = clusters[i]
        if len(spiketimes) == 0:
            print('No spikes in this cluster')
        else:
            storespikes['base_0'] = spiketimes[np.where((spiketimes >= stimtp[0][0]-60*Fs)& (spiketimes < stimtp[0][0]) )] - (stimtp[0][0]-60*Fs)
            storespikes['dur_0'] = spiketimes[np.where((spiketimes >= (stimtp[0][0] )) & (spiketimes < (stimtp[0][1])))[0]] - (stimtp[0][0])
            storespikes['post_0'] = spiketimes[np.where((spiketimes >= (stimtp[0][1])) & (spiketimes < stimtp[0][1] + 60*Fs))[0]] - (stimtp[0][1])
            
            for stimint in range(len(stimtp) -1):
                storespikes['base_' + str(stimint+1)] = spiketimes[np.where((spiketimes >= stimtp[stimint+1][0] + segmentlengths[0][stimint] -60*Fs ) & (spiketimes < (stimtp[stimint+1][0] + segmentlengths[0][stimint] )))[0]] - (stimtp[stimint+1][0] + segmentlengths[0][stimint] -60*Fs)
                storespikes['dur_' + str(stimint+1)] = spiketimes[np.where((spiketimes >= (stimtp[stimint+1][0] + segmentlengths[0][stimint] )) & (spiketimes < (stimtp[stimint+1][1] + segmentlengths[0][stimint] )))[0]] - (stimtp[stimint+1][0] + segmentlengths[0][stimint])
                storespikes['post_' + str(stimint+1)] = spiketimes[np.where((spiketimes >= (stimtp[stimint+1][1] + segmentlengths[0][stimint] )) & (spiketimes < (stimtp[stimint+1][1] + segmentlengths[0][stimint] +60*Fs)))[0]] - (stimtp[stimint+1][1] + segmentlengths[0][stimint])


            
            
            numprewindows = np.zeros(len(stimtp))
            numdurwindows = np.zeros(len(stimtp))
            
            
            
            for b in range(len(stimtp)):
                numdurwindows[b] = (((stimtp[b][1]) - (stimtp[b][0]))/Fs)/10
                numprewindows[b] = (stimtp[b][0]/Fs)/10
            
            totwindows = np.add(numdurwindows,numprewindows)
            
            preFRstack = []
            
            
            totstack = []
            durstack = []
            prestack = []
            for b in range(len(stimtp)):
                tempdurstack = []
                tempprestack = []
                for count in range(int(np.floor(numdurwindows[b]))):
                    temparray = storespikes['dur_' + str(b)][np.where((storespikes['dur_'+str(b)] > (count*10*Fs)) & (storespikes['dur_'+str(b)] < ((count+1)*10*Fs) ))[0]]
                    tempdurstack.append(len(temparray)/10)
                for count in range(int(np.floor(numprewindows[b]))):
                    temparray = storespikes['base_' + str(b)][np.where((storespikes['base_'+str(b)] > (count*10*Fs)) & (storespikes['base_'+str(b)] < ((count+1)*10*Fs) ))[0]]
                    tempprestack.append(len(temparray)/10)
                tempdurstack = np.array(tempdurstack)
                tempprestack = np.array(tempprestack)
                durstack.append(np.array(tempdurstack))
                prestack.append(np.array(tempprestack))
                FRstack = np.concatenate((tempprestack,tempdurstack))
                totstack.append(FRstack)
            durstack = np.array(durstack)    
            prestack = np.array(prestack)    
            

                    
            for b in range(len(stimtp)):
                if (len(np.where(totstack[b] >= 0.25)[0]) > totwindows[b]/2) & (len(np.where(prestack[b] >= 0.25)[0]) >= 1) & (len(np.where(durstack[b] >= 0.25)[0]) >= 1):
                    print('Timepoint:' +str(b))
                    waveform = (stimwaveform[b][0] - stimwaveform[b][1])[stimtp[b][0]:stimtp[b][1]]
                    [basephasespike, baseampspike] = phasescalc(waveform, storespikes['base_' + str(b)])
                    
                    [durphasespike, durampspike] = phasescalc(waveform, storespikes['dur_' + str(b)])
                    
                    [postphasespike, postampspike] = phasescalc(waveform, storespikes['post_' + str(b)])
                    
                    allphases.append([i, b, basephasespike, durphasespike, postphasespike, amp[b], freq[b], storespikes['base_' + str(b)], storespikes['dur_' + str(b)], storespikes['post_' + str(b)], waveform, clusterlabel])
                    tpcounter[b] += 1            
                    #print('Done')
                else:
                    #print('Not enough spikes')
                    continue
                
                
        
    with open(folder +"1minsACS_allphases", "wb") as fp:   #Pickling
        pickle.dump([allphases, len(clusters), tpcounter], fp)
        

#%%

#Insert location of 20-min sACS experiment folders   
directories = sorted(glob.glob('/20minsACSdata/202*/'))
for directory in directories: 
    print(directory)
    stimdata = np.load(directory + 'stimtimepoints.npy')
    MUAdata = mat73.loadmat(sorted(glob.glob(directory + 'sACS*/**/'))[0] + 'MUAfiltdata.mat')
    segmentlengths = np.load(directory + 'segmentlengths.npy')

    amp = MUAdata['MUAsave']['amp']
    freq = MUAdata['MUAsave']['freq']
    Fs = MUAdata['MUAsave']['Fs']
    stimon = stimdata[0]*Fs
    stimoff = stimdata[1]*Fs
    stimwaveform = (MUAdata['MUAsave']['stimdata'][0] - MUAdata['MUAsave']['stimdata'][1])[int(stimon):int(stimoff)]
    
    
    clusterfile = directory + 'MUAfiltdata/MUAfiltdata.result-merged.hdf5'
    c = h5py.File(clusterfile, 'r')['spiketimes'] 
    clusters = list(c.keys())


    
    
    allphases = []
    allphases15 = []
    
    for i in range(len(segmentlengths)-1):
        segmentlengths[i+1]= segmentlengths[i+1] + segmentlengths[i]
    
    print('post stim recording is ' +str(((segmentlengths[0] - stimoff)/Fs)/60) + ' min long')
    
    if stimon > 20*60*Fs:
        start2clip = stimon-(20*60*Fs)
    else:
        start2clip = 0

    i = []
    for i in range(len(clusters)):
        spiketimes = c[clusters[i]]
        if len(spiketimes) == 0:
            print('No spikes in this cluster')
        else:
            storespikes = {}
            storespikes['base'] = spiketimes[np.where((spiketimes >= start2clip) &(spiketimes < stimon))[0]] - start2clip
            storespikes['dur'] = spiketimes[np.where((spiketimes  > (stimon)) & (spiketimes < stimoff))[0]] - stimon
            storespikes['post_0'] = spiketimes[np.where((spiketimes > stimoff)&(spiketimes < segmentlengths[0]))[0]] - stimoff
            
            for seg in range(len(segmentlengths) - 1):
                storespikes['post_' + str(seg+1)] =  spiketimes[np.where((spiketimes > segmentlengths[seg])&(spiketimes < segmentlengths[seg+1]))[0]] - segmentlengths[seg]
                
                
            numprewindows = int(np.floor(((stimon - start2clip)/Fs)/10))
            numdurwindows = int(np.floor(((stimoff-stimon)/Fs)/10))
            totwindows = numprewindows + numdurwindows
            
            
            preFRstack = []
            durFRstack = []
            for count in range(numprewindows):
                temparray = storespikes['base'][np.where((storespikes['base'] > (count*10*Fs))&(storespikes['base'] < ((count+1)*10*Fs)))]
                preFRstack.append(len(temparray)/10)
                
            preFRstack = np.array(preFRstack)
            
            
            for count in range(numdurwindows):
                temparray = storespikes['dur'][np.where((storespikes['dur'] > (count*10*Fs))&(storespikes['dur']< ((count+1)*10*Fs)))]
                durFRstack.append(len(temparray)/10)
                
            durFRstack = np.array(durFRstack)
            
            FRstack = np.concatenate((preFRstack,durFRstack))
            

        
        #if ((len(prestim)/(stimon/Fs) > 0.25) and (len(durstim)/((stimoff - (stimon+2*Fs))/Fs)>0.25)):
            if len(np.where(FRstack >= 0.25)[0]) > (totwindows/2) and len(storespikes['base'])>2 and len(storespikes['dur'])>2:   
                
                
                [phasespike, ampspike] = phasescalc(stimwaveform, storespikes['base'])
        #Might need to change this - check what info is needed downstream
                allphases.append([phasespike,i,0,amp,freq,(stimon - start2clip)/Fs,storespikes['base'], clusters[i]])
        
                [phasespike, ampspike] = phasescalc(stimwaveform,storespikes['dur'])
                allphases.append([phasespike,i,1,amp,freq,(stimoff-stimon)/Fs,storespikes['dur'], clusters[i]])
                
                [phasespike, ampspike] = phasescalc(stimwaveform,storespikes['post_0'])
                allphases.append([phasespike,i,2,amp,freq,(segmentlengths[0] - stimoff)/Fs,storespikes['post_0'], clusters[i]])
                
                for seg in range(len(segmentlengths)-1):
                    [phasespike, ampspike] = phasescalc(stimwaveform,storespikes['post_' + str(seg+1)])
                    allphases.append([phasespike,i,2+seg,amp,freq,(segmentlengths[seg+1] - segmentlengths[seg])/Fs,storespikes['post_' + str(seg+1)], clusters[i]])
        
            
            else:
                print('Not enough spikes')
        
        
    with open(directory + "20minsACS_allphases", "wb") as fp:   #Pickling
        pickle.dump(allphases, fp)

        
        
        
        
        
        
