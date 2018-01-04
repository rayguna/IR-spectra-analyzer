# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:17:39 2018

@author: rayg

Stores common functions:
    1. standardize units and normalize 
    2. Interpolate spectrum
    
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import my_jcamp as jcamp


#1. standardize units and normalize 
def standardize_units_and_normalize_peaks(jcamp_dict):
        """Given a jcamp dictionary, standardize x and y units to 1/cm and absorbance
           , respectively and return a modified dictionary. 
           Args:
           jcamp_dict=a jcamp dictionary    
           Returns:
           jcamp_dict=a modified jcamp dictionary   
        """

        if jcamp_dict['yunits'] == "ABSORBANCE":
            pass 

        elif jcamp_dict['yunits'] == "TRANSMISSION" or "TRANSMITTANCE":
            jcamp_dict['y'] = 2 - np.log10(100*(jcamp_dict['y']+0.1)) #avoid division by zero
            #normalize
            jcamp_dict['y'] = (jcamp_dict['y']-min(jcamp_dict['y'])) / max(jcamp_dict['y'])



        # check xunits: if in microns, change to 1/cm
        if jcamp_dict['xunits'] == "MICROMETERS":
            jcamp_dict['x'] = 10000 / jcamp_dict['x']

        # uniformize data, #2:
        # normalize absorbance peaks (y-values) to between 0 and 1.
        jcamp_dict['y'] = (jcamp_dict['y']-min(jcamp_dict['y'])) / max(jcamp_dict['y'])

        return jcamp_dict

#2. Interpolate spectrum
def interpolate_spectrum(jcamp_dict,xmin=800,xmax=3000,res=0.5):
    """Given a jcamp dictionary, interpolate the data points to the specified 
       xmin, xmas, and resolution.
       Args:
       jcamp_dict=a jcamp dictionary    
       Returns:
       xmin=minimum x-value
       xmax=maximum x-value
       res=x-inteval    
    """
    
    #standardize units
    my_dict=standardize_units_and_normalize_peaks(jcamp_dict)
    
    #interpolate data set 1
    dx=(my_dict['lastx']-my_dict['firstx'])/my_dict['npoints']
    
    my_list_x=[] #generate x values
    for i in range(my_dict['npoints']):
        my_list_x.append(my_dict['firstx']+(dx*i))
        
    tck = interpolate.splrep(my_list_x, my_dict['y'], s=0) # s is smoothing
    
    jcamp_dict['x'] = np.arange(xmin,xmax,res) #the second limit is excluded.
    jcamp_dict['y'] = interpolate.splev(jcamp_dict['x'], tck, der=0) #der refers to the derivative
    
    return jcamp_dict

#3. Calculate Euclidean distance
def calculate_euclidean_dist(spectrum1,spectrum2):
    """Given two spectra, calculate the Euclidean distance
       Args:
       jcamp1=spectrum 1 as a jcamp dictionary
       jcamp2=spectrum 2 as a jcamp dictionary
       Returns:
       dst=the Euclidean distance 
       a plot showing the two spectra and their difference
    """     
    
    #standardize units
    spectrum1=standardize_units_and_normalize_peaks(spectrum1)
    spectrum2=standardize_units_and_normalize_peaks(spectrum2)
    
    #determine the xlimits to use and interpolate data accordingly
    first_x=max([spectrum1['x'][0],spectrum2['x'][0]])
    last_x=min([spectrum1['x'][-1],spectrum2['x'][-1]])
    
    #interpolate data set 
    spectrum1=interpolate_spectrum(spectrum1,first_x,last_x)
    spectrum2=interpolate_spectrum(spectrum2,first_x,last_x)
    
    #determine the xlimits to use and interpolate data accordingly
    plt.plot(spectrum1['x'],spectrum1['y'],'b',label=spectrum1['title'])
    plt.plot(spectrum2['x'],spectrum2['y'],'r',label=spectrum2['title'])
    
    #subtract y values of the two plots (i.e., calculate the difference between the two plots)
    y_diff=spectrum1['y']-spectrum2['y']
    #plot the difference
    plt.plot(spectrum1['x'],y_diff,'g--', label='difference')
    
    plt.legend()
    plt.xlabel('Wavenumber (cm-1)')
    plt.ylabel('Absorbance')
    plt.show()
    
    ##########################################################################################
    #Calculate correlation coefficient
    #print(np.corrcoef(spectrum1['y'],spectrum2['y'])[0,1])
    
    #calculate euclidean distance
    #a. using scipy
    from scipy.spatial import distance
    dst1 = distance.euclidean(spectrum1['y'],spectrum2['y'])
    ##print(dst1)
    
    #b. using numpy
    #dst2 = np.linalg.norm(spectrum1['y']-spectrum2['y'])
    
    return dst1
    #print('Euclidean distance: %.2f' %(dst2))

"""
#TEST    
#read jcamp of two spectra data sets
spectrum1=jcamp.JCAMP_reader("1-Ethoxy-4-nitrobenzene_100-29-8") #1st data set
spectrum2=jcamp.JCAMP_reader("Acetanilide 2-chloro-4-tert-butyl-_100141-30-8") #2nd data set

calculate_euclidean_dist(spectrum1,spectrum2)

"""



    


