# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 12:27:30 2023

@author: Anang
"""

'''
This is master python file to do the image analysis with stepwise. Follow the steps 
'''
import sys
from termcolor import colored
sys.exit()

# exec(open('./scripts/seg_phase_imges.py').read())



#change the folder name in vars_allround1data

#step0: create the folders for output analysis
#%%

exec(open('./scripts/_runFirst.py').read())

#%%
#step1: segment the bright field cell images using Yeaz

exec(open('./scripts/segment_cell.py').read())
#%%
#step2: Unmixed the organlles or pre-processes

'''
peroxisome and vacule 
'''
exec(open('./scripts/preprocess_blueFinal.py').read())

'''
ER 
'''
exec(open('./scripts/preprocess_green.py').read())
'''
golgi; mitochondria; LD
'''
exec(open('./scripts/preprocess_yellowNredFinal.py').read())

print("======================================================= ")
print("====",colored("Go for ilastik segmentation","magenta",attrs=["bold"]),"===========")
print("======================================================= ")

# quit()
sys.exit()

###========================================================================####
###========================================================================####

#Step3: Go for ilastik segmentation
'''
use the folder name: ilastk_control_libfile
'''
###========================================================================####
###========================================================================####


#Step4: Check the spectra of red and blue color organelles



exec(open('./scripts/raw2spectra_blue.py').read())

#%%
exec(open('./scripts/raw2spectra_red.py').read())

#%%
# exec(open('./scripts/raw2spectra.py').read()) #for all folder together 


#step5: label the images or post-process

'''
    "peroxisome",
"LD"
"golgi"

'''
exec(open('./scripts/postprocess_globular.py').read())

'''
ER
'''
exec(open('./scripts/postprocess_ER.py').read())

'''
Mitochondrria
'''
exec(open('./scripts/postprocess_mito.py').read())
'''
Vacuoles
'''
exec(open('./scripts/postprocess_vacuole.py').read())

#step6: Go for measurement

exec(open('./scripts/measure_cell.py').read())

exec(open('./scripts/measure_organelle.py').read())

# exec(open('./scripts/measure_organelle_testIntensity.py').read())

print("======================================================= ")
print("====",colored("Go for analysis of the csv file generated","magenta",attrs=["bold"]),"===========")
print("======================================================= ")


####=================================================================
print("====================",colored("read files: use csv2figure_new.py","red"),"===================================")
####=================================================================

