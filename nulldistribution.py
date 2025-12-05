#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##
##After extracting the phases of spikes times corresponding to the stimulation waveform, this code
##shuffles the phases 10 000 times in baseline and during/post stimulation to create a null distribution
## This null distribution will then be used against actual differences in phases/# of spikes to assess 
## which neurons were significantly affected by sACS

#%%

import pandas as pd
import numpy as np
import math
import pickle
import glob
from multiprocessing import Pool
import random


#%%
#Insert location of 1-min sACS experiment folders
files = glob.glob('/1minsACS_allphases')

n = 10000 #Number of times to repeat iteration
Fs = 24414.0625
binsize = int(0.1*Fs)  # 100ms timebins
clip = 5 #Clipping first 5 seconds after stimulation onset to avoid errors arond edges of artefact removal

numbins = int(np.floor((55*Fs)/binsize))

def PLVcalc(phases):
    if len(phases) >1:
        PLVval = abs(sum(math.e**(1j*phases))/len(phases))
    else:
        PLVval = np.nan
    return PLVval

#Have to check this one as I'm not sure if it will work
def PPCcalc(phasespike):
    phasespike = np.float128(phasespike) 
    # ip = argparse.ArgumentParser()
    # ip.add_argument('phasespike')
    N = len(phasespike)
    if N > 1:
        scale_factor = (2/(N*(N-1)))
        sphases = np.sin(phasespike)
        cphases = np.cos(phasespike)
    
        ppc0 = np.zeros(N-1)
        ii = []
        for ii in range(0,N-1):
            ppc0[ii] = np.sum(((cphases[ii]*cphases[(ii+1):N]) + (sphases[ii]*sphases[(ii+1):N])))
            # def innerppc(ii):
                #     ppctemp = np.sum(((cphases[ii]*cphases[(ii+1):N]) + (sphases[ii]*sphases[(ii+1):N])))
                #     return ppctemp
                # for ii in pool.imap(innerppc,range(0,N-1)):
                    #     ppc0[ii] = ii
        ppc0 = np.sum(ppc0) * scale_factor
    else:
        ppc0 = float('nan')
    return ppc0

def min1_calcnulldist(neur):
    nullPLV = []
    nullPPC = []
    nullPLVdelt = []
    nullPPCdelt = []
    nullFRdelt = []
    
    neurdata = files[neur]
    
    combine1 = neurdata[2][np.where(neurdata[7] > clip*Fs)]
    combine2 = neurdata[3][np.where(neurdata[8] > clip*Fs)]
    
        
    tempPLV = []
    tempPPC = []
    i = []
    
    prebins = []
    durbins = []
    totbins = []
    for i in range(numbins):
        prebins.append(np.size(
            np.where((neurdata[7] > (i*binsize + 5*Fs))&(neurdata[7] <= (i*binsize + binsize + 5*Fs)))[0]))
        durbins.append(np.size(
            np.where((neurdata[8] > (i*binsize + 5*Fs))&(neurdata[8] <= (i*binsize + binsize + 5*Fs)))[0]))
        
    totbins = np.concatenate((prebins, durbins))
    
    for i in range(n):
        if len(combine1) > len(combine2):
            random.shuffle(combine1)
            phasesshuf = np.concatenate((combine1[0:len(combine2)], combine2))
        elif len(combine1) <= len(combine2):
            random.shuffle(combine2)
            phasesshuf = np.concatenate((combine1, combine2[0:len(combine1)]))
        random.shuffle(phasesshuf)
        temppre = phasesshuf[0:int(len(phasesshuf)/2)]
        tempdur = phasesshuf[int(len(phasesshuf)/2):]
        tempPLV.append([PLVcalc(temppre),PLVcalc(tempdur)])
        tempPPC.append([PPCcalc(temppre),PPCcalc(tempdur)])
        random.shuffle(totbins)
        temppre = np.mean(totbins[0:int(len(totbins)/2)])/0.1
        tempdur = np.mean(totbins[int(len(totbins)/2):])/0.1
        nullFRdelt.append(tempdur-temppre)
    nullPLV = np.array(tempPLV)
    nullPPC = np.sqrt(abs(np.array(tempPPC)))
    nullPLVdelt.append([nullPLV[:,1] - nullPLV[:,0]])
    nullPPCdelt.append([nullPPC[:,1] - nullPPC[:,0]])
    
    return [nullPLVdelt, nullPPCdelt, nullFRdelt, neur, neurdata[0], neurdata[1], neurdata[5], neurdata[6]]
    


def min20_calcnulldist_dur(nulldataval):
       nullPLV = []
       nullPPC = []
       nullPLVdelt = []
       nullPPCdelt = []
       nullFRdelt = []
       lenbaseline = nullphasedata[int(2*nulldataval)][5]
       
       
       combine1 = nullphasedata[int(2*nulldataval)][0][np.where((nullphasedata[int(2*nulldataval)][6] >= (lenbaseline*Fs - 20*60*Fs + clip*Fs)) & (nullphasedata[int(2*nulldataval)][6] < (lenbaseline*Fs)))[0]]
       combine2 = nullphasedata[int(2*nulldataval +1)][0][np.where((nullphasedata[int(2*nulldataval+1)][6] >clip*Fs)&(nullphasedata[int(2*nulldataval +1)][6] < 20*60*Fs))[0]]
       
  
       
       
       tempPLV = []
       tempPPC = []
       i = []
       
       prebins = []
       durbins = []
       totbins = []
       for i in range(numbins):
           prebins.append(np.size(
               np.where((nullphasedata[int(2*nulldataval)][6] > (i*binsize + clip*Fs))&(nullphasedata[int(2*nulldataval)][6] <= (i*binsize + binsize + clip*Fs)))[0]))
           durbins.append(np.size(
               np.where((nullphasedata[int(2*nulldataval + 1)][6] > (i*binsize + clip*Fs)) & (nullphasedata[int(2*nulldataval+1)][6] <= (i*binsize + binsize + clip*Fs)))[0]))
       totbins = np.concatenate((prebins, durbins))
       
       
       for i in range(n):
           if len(combine1) > len(combine2):
                random.shuffle(combine1)
                phasesshuf = np.concatenate((combine1[0:len(combine2)], combine2))
           elif len(combine1) <= len(combine2):
                random.shuffle(combine2)
                phasesshuf = np.concatenate((combine1, combine2[0:len(combine1)]))
           random.shuffle(phasesshuf)
           temppre = phasesshuf[0:int(len(phasesshuf)/2)]
           tempdur = phasesshuf[int(len(phasesshuf)/2):]
           tempPLV.append([PLVcalc(temppre),PLVcalc(tempdur)])
           tempPPC.append([PPCcalc(temppre),PPCcalc(tempdur)])
           random.shuffle(totbins)
           temppre = np.nanmean(totbins[0:int(len(totbins)/2)])/0.1
           tempdur = np.nanmean(totbins[int(len(totbins)/2):])/0.1
           nullFRdelt.append(tempdur-temppre)
       nullPLV = np.array(tempPLV)
       nullPPC = np.sqrt(abs(np.array(tempPPC)))
       nullPLVdelt.append([nullPLV[:,1] - nullPLV[:,0]])
       nullPPCdelt.append([nullPPC[:,1] - nullPPC[:,0]])
       
       return [nullPLVdelt, nullPPCdelt, nullFRdelt, nulldataval]


def min20_calcnulldist_post(baselocation):
    nullPLV = []
    nullPPC = []
    nullPLVdelt = []
    nullPPCdelt = []
    nullFRdelt = []
       
    lenbaseline = phasedata[baselocation][5]
    lenpoststim1 = phasedata[baselocation+2][5]
    

    checkimprestim = 0
    
    if lenpoststim1 >= timewin:
        combine2 = phasedata[baselocation+2][0][np.where(phasedata[baselocation+2][6] <= timewin*Fs)[0]]
        posttp = phasedata[baselocation+2][6][np.where(phasedata[baselocation+2][6] <= timewin*Fs)[0]]
        
        
        if lenbaseline >= timewin + 60:
            combine1 = phasedata[baselocation][0][np.where((phasedata[baselocation][6] > 60*Fs)&(phasedata[baselocation][6] <=  (60*Fs + timewin*Fs)))[0]]
            combine11 = phasedata[baselocation][0][np.where((phasedata[baselocation][6] >= (lenbaseline*Fs) - timewin*Fs)&(phasedata[baselocation][6] <  (lenbaseline*Fs)))[0]]
            basetp = phasedata[baselocation][6][np.where((phasedata[baselocation][6] > 60*Fs)&(phasedata[baselocation][6] <=  (60*Fs + timewin*Fs)))[0]]
            basetp11 =  phasedata[baselocation][6][np.where((phasedata[baselocation][6] >= (lenbaseline*Fs) - timewin*Fs)&(phasedata[baselocation][6] <  (lenbaseline*Fs)))[0]]

            checkimprestim = 1
            
            
        else:
            combine1 = phasedata[baselocation][0][np.where(phasedata[baselocation][6] <= timewin*Fs)[0]]
            basetp = phasedata[baselocation][6][np.where(phasedata[baselocation][6] <= timewin*Fs)[0]]
        
        numbins = int((timewin*Fs)/binsize)
        
        
        
    elif lenpoststim1 < timewin:
        print('First post stim period is only: ' + str(lenpoststim1))
        print('Resetting timewin to match')
        
        timewin2 = lenpoststim1
        
        combine2 = phasedata[baselocation+2][0][np.where(phasedata[baselocation+2][6] <= timewin2*Fs)[0]]
        posttp = phasedata[baselocation+2][6][np.where(phasedata[baselocation+2][6] <= timewin2*Fs)[0]]
        
        
        
        if lenbaseline >= (timewin2+60):
            combine1 = phasedata[baselocation][0][np.where((phasedata[baselocation][6] > 60*Fs)&(phasedata[baselocation][6] <=  (60*Fs + timewin2*Fs)))[0]]
            combine11 = phasedata[baselocation][0][np.where((phasedata[baselocation][6] >= (lenbaseline*Fs) - timewin2*Fs)&(phasedata[baselocation][6] <  (lenbaseline*Fs)))[0]]
            basetp = phasedata[baselocation][6][np.where((phasedata[baselocation][6] > 60*Fs)&(phasedata[baselocation][6] <=  (60*Fs + timewin2*Fs)))[0]]
            basetp11 =  phasedata[baselocation][6][np.where((phasedata[baselocation][6] >= (lenbaseline*Fs) - timewin2*Fs)&(phasedata[baselocation][6] <  (lenbaseline*Fs)))[0]]

            checkimprestim = 1

        else:
            combine1 = phasedata[baselocation][0][np.where(phasedata[baselocation][6] <= timewin2*Fs)[0]]    
            basetp = phasedata[baselocation][6][np.where(phasedata[baselocation][6] <= timewin2*Fs)[0]]
        
        numbins = int(np.floor(timewin2*Fs)/binsize)
    
    if checkimprestim == 0:
                 
        if len(combine1) > len(combine2): 
            phasesshuf = np.concatenate((combine1[0:len(combine2)], combine2))
        elif len(combine1) <= len(combine2): 
            phasesshuf = np.concatenate((combine1, combine2[0:len(combine1)]))
        tempPLV = []
        tempPPC = []
        i = []
           
        prebins = []
        durbins = []
        totbins = []
        for i in range(numbins):
            prebins.append(np.size(
                   np.where((basetp > ( i*binsize))&(basetp <= ( i*binsize + binsize)))[0]))
            durbins.append(np.size(
                   np.where((posttp > (i*binsize)) & (posttp <= (i*binsize + binsize)))[0]))
        totbins = np.concatenate((prebins, durbins))
            
            
           
           
        for i in range(n):
            random.shuffle(phasesshuf)
            temppre = phasesshuf[0:int(len(phasesshuf)/2)]
            tempdur = phasesshuf[int(len(phasesshuf)/2):]
            tempPLV.append([PLVcalc(temppre),PLVcalc(tempdur)])
            tempPPC.append([PPCcalc(temppre),PPCcalc(tempdur)])
            random.shuffle(totbins)
            temppre = np.nanmean(totbins[0:int(len(totbins)/2)])/0.1
            tempdur = np.nanmean(totbins[int(len(totbins)/2):])/0.1
            nullFRdelt.append(tempdur-temppre)
        nullPLV = np.array(tempPLV)
        nullPPC = np.sqrt(abs(np.array(tempPPC)))
        nullPLVdelt.append([nullPLV[:,1] - nullPLV[:,0]])
        nullPPCdelt.append([nullPPC[:,1] - nullPPC[:,0]])

    elif checkimprestim == 1:
        if len(combine11) > len(combine2): 
            phasesshuf = np.concatenate((combine11[0:len(combine2)], combine2))
        elif len(combine11) <= len(combine2): 
            phasesshuf = np.concatenate((combine11, combine2[0:len(combine11)]))
        tempPLV = []
        tempPPC = []
        i = []
           
        prebins = []
        durbins = []
        totbins = []
        for i in range(numbins):
            prebins.append(np.size(
                   np.where((basetp11 > (i*binsize))&(basetp11 <= (i*binsize + binsize)))[0]))
            durbins.append(np.size(
                   np.where((posttp > (i*binsize)) & (posttp <= (i*binsize + binsize)))[0]))
        totbins = np.concatenate((prebins, durbins))
            
            
           
           
        for i in range(n):
            random.shuffle(phasesshuf)
            temppre = phasesshuf[0:int(len(phasesshuf)/2)]
            tempdur = phasesshuf[int(len(phasesshuf)/2):]
            tempPLV.append([PLVcalc(temppre),PLVcalc(tempdur)])
            tempPPC.append([PPCcalc(temppre),PPCcalc(tempdur)])
            random.shuffle(totbins)
            temppre = np.mean(totbins[0:int(len(totbins)/2)])/0.1
            tempdur = np.mean(totbins[int(len(totbins)/2):])/0.1
            nullFRdelt.append(tempdur-temppre)
        nullPLV = np.array(tempPLV)
        nullPPC = np.sqrt(abs(np.array(tempPPC)))
        nullPLVdelt.append([nullPLV[:,1] - nullPLV[:,0]])
        nullPPCdelt.append([nullPPC[:,1] - nullPPC[:,0]])
        
    return [nullPLVdelt, nullPPCdelt, nullFRdelt, baselocation]
    
    
    


#%%
for file in files:
    beep = pd.read_pickle(file)
    files = beep[0] #phase data for all neurons
    
    
    nullfinal = []

    if __name__ == '__main__':  
        p = Pool()
        nullfinal.append(p.map(min1_calcnulldist,range((len(files)))))
        
    with open(file[0:file.find('/1minsACS') +1] + "1minsACS_nulldistribution", "wb") as fp:   #Pickling
        pickle.dump([nullfinal], fp)


#%%
#Computing null distribution for 20-min sACS experiments during 20-minutes of stimulation delivery
#Insert location of 20-min sACS experiment folders
files = glob.glob('/20minsACS_allphases')


for file in files:
    print(file)
    phasedata = pd.read_pickle(file)
    

    #Collect the pre and during stimulation periods and combine into a new dataset
    nullphasedata = []
    i = []
    for i in range(len(phasedata)):
        if (phasedata[i][2] == 0) or (phasedata[i][2] == 1):
            nullphasedata.append(phasedata[i])
        else:
            continue

    nullfinal = []

    if __name__ == '__main__':  
        p = Pool()
        nullfinal.append(p.map(min20_calcnulldist_dur,range(int(len(nullphasedata)/2))))


    with open(file[0:file.find('202')+11] + "20minsACS_nulldistribution_durstim", "wb") as fp:   #Pickling
        pickle.dump([nullfinal], fp)




#%%
#Computing null distribution for 20-min sACS experiments in the 10minutes following stimulation offset
timewin = 60*10


for file in files:
    print(file)
    phasedata = pd.read_pickle(file)


    #Collect the pre and during stimulation periods and combine into a new dataset
    baseloc = []
    i = []
    for i in range(len(phasedata)):
        if (phasedata[i][2] == 0) :
            baseloc.append(i)
        else:
            continue
    baseloc = np.array(baseloc)

    nullfinal = []
    
    
    
    if __name__ == '__main__':  
        p = Pool()
        nullfinal.append(p.map(min20_calcnulldist_post,baseloc))


    with open(file[0:file.find('202')+11] + "20minsACS_nulldistribution_poststim", "wb") as fp:   #Pickling
        pickle.dump([nullfinal], fp)
    
        
        
        












