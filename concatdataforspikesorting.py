#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:54:09 2022

@author: mbellvila
"""

#After running MUAartrem.m to clean the stimulation artefact from the MUA data, 
#this script concatenates data from rat subfolders to subsequently analyze using SpyKING Circus

#%%
import glob
import numpy as np
import scipy.io
import mat73
import tdt

#%%
#This should be edited to the folder containing all 1-min sACS experiments
#This section of the code will concatenate and save all the filtered MUA data for each rat of the 1-min sACS experiments
dirstosearch = glob.glob('/locationof1minACSdata/202*/')

for directory in dirstosearch:
    subfolders = sorted(glob.glob(directory + 'sACS*/**/*MUAfiltdata.mat'))
    
    combinedMUA = []
    temp = mat73.loadmat(subfolders[0])['MUAsave']
    combinedMUA = temp['MUAdata']
    lengthofsegment = [len(combinedMUA[0])]
    stimtp = [temp['stimon']]
    stimwf = [temp['stimdata']]
    stimfreq = [temp['freq']]
    stimamp = [temp['amp']]
    for i in range(len(subfolders)-1):
        temp = mat73.loadmat(subfolders[i+1])['MUAsave']
        data = temp['MUAdata']
        stimtp.append(temp['stimon'])
        combinedMUA = np.concatenate((combinedMUA,data), axis =1)
        lengthofsegment.append(len(data[0]))
        stimwf.append(temp['stimdata'])
        stimfreq.append(temp['freq'])
        stimamp.append(temp['amp'])
        
        
    np.save(directory + 'MUAfiltdata.npy',combinedMUA)    
    np.save(directory + 'segmentlengths.npy',np.array([lengthofsegment, stimtp, stimwf, stimfreq, stimamp], dtype=object))
    
    
#%%
#This should be edited to the folder containing all 20-min sACS experiments
#This concatenates all the MUA data per rat.
directories = glob.glob('/locationof20minACSdata/202*/')
for directory in directories:
    print(directory)
    filestocombine = sorted(glob.glob(directory + 'sACS*/**/*MUAfiltdata.mat'))
    combinedMUA = mat73.loadmat(filestocombine[0])['MUAsave']['MUAdata']
    lengthofsegment = [len(combinedMUA[0])]
    for i in range(len(filestocombine)-1):
        temp = mat73.loadmat(filestocombine[i+1])['MUAsave']
        data = temp['MUAdata']
        combinedMUA = np.concatenate((combinedMUA,data),axis=1)
        lengthofsegment.append(len(data[0]))
    np.save(directory + 'MUAfiltdata.npy', combinedMUA)
    np.save(directory + 'segmentlengths.npy',lengthofsegment)    

#This portion of the code saves the stimulation timepoints for 20-min sACS experiments
for directory in directories:
    print(directory)
    data = tdt.read_block(sorted(glob.glob(directory + 'sACS*/**/'))[0])
    
    stimon = data.epocs.Wfrq.onset[0]
    stimoff = data.epocs.Wfrq.offset[0]

    print(str(stimon))
    print(str(stimoff))

    np.save(directory + 'stimtimepoints.npy', [stimon, stimoff] )



















