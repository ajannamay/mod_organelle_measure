# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:57:56 2023

@author: Anang
"""

import organelle_measure.vars_allround1data as var
import os



list_folders = var.list_folders
isExist = os.path.exists(f'{var.path_out}')
if not isExist:
    print("exist")
    os.mkdir(f'{var.path_out}')
for folders in list_folders:
    os.mkdir(f'{var.path_out}/{folders}')

for folders in list_folders:
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["segmented_cell"]}')
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["organelle_unmixed"]}')
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["ilastk_segmentated"]}')
    
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["label_organelle"]}')
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["spectra"]}')
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["measurements"]}')
    os.mkdir(f'{var.path_out}/{folders}/{var.output_folders["ilastk_control_libfile"]}') 
