# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:54:53 2023

@author: Anang
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:37:21 2023

@author: Anang
"""

# This script prints out the spectrum of the raw images,
# in order to validate our unmixing.

# This file assumes running on the root dir of this project!

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import xmltodict
from pathlib import Path
from copy import deepcopy
from skimage import io,util
from organelle_measure.tools import load_nd2_plane,get_nd2_size
import organelle_measure.vars_allround1data as var

def read_spectra_xml(filepath):
    with open(filepath,'rb') as file_xml:
        data_dict = xmltodict.parse(file_xml.read())
    # get the effective key under 'variant'
    count_variant = 0
    for key_variant in data_dict['variant'].keys():
        count_variant += 1
    assert count_variant==2, "Wrong Length of Variant's Children!"

    # get the keys of each individual spectrum
    keys_unmix = []
    count_unmix = 0
    for key_unmix in data_dict["variant"][key_variant].keys():
        if key_unmix[0] != '@':
            count_unmix += 1
            keys_unmix.append(key_unmix)

    # get the data points
    num_points = int(data_dict["variant"][key_variant][key_unmix]["Spectrum"]["uiCount"]["@value"])
    properties = ['eType', 'uiWavelength', 'dWavelength','dTValue']
    node = data_dict["variant"][key_variant][key_unmix]["Spectrum"]["pPoint"][f"Point0"]
    dict_type = {prop:int if "int" in node[prop]["@runtype"] else float for prop in properties}
    list_spectra = []
    for key_unmix in keys_unmix:
        list_spectra.append({prop:np.zeros((num_points,),dtype=dict_type[prop]) for prop in properties})
        for point in range(num_points):
            node = data_dict["variant"][key_variant][key_unmix]["Spectrum"]["pPoint"][f"Point{point}"]
            for key_spectra in list_spectra[-1].keys():
                type_point = int if "int" in node[key_spectra]["@runtype"] else float # how can we not avoid this
                list_spectra[-1][key_spectra][point] = type_point(node[key_spectra]["@value"])
    for i,spectrum in enumerate(list_spectra):
        spectrum["unmix"] = i
    df = pd.concat((pd.DataFrame(spectrum) for spectrum in list_spectra),ignore_index=True)
    return df

def get_spectra_img(file1,file2,file_raw):
    img_label1 = io.imread(str(file1)) #ilatic seg images organeel first
    img_label2 = io.imread(str(file2))#second orgennels
    img_raw    = load_nd2_plane(str(file_raw),frame="czyx",axes="t",idx=0) #two orgenells together
    # shape = get_nd2_size(str(file_raw))

    label1 = (img_label1>1)
    label2 = (img_label2>1)
    # label0 = np.invert(label1)*np.invert(label2)

    # spectra0 = np.array(img_raw[:,label0]).transpose()
    spectra1 = np.array(img_raw[:,label1]).transpose()
    spectra2 = np.array(img_raw[:,label2]).transpose()

    # max0 = np.max(np.array(spectra0),axis=1)[:,np.newaxis]
    max1 = np.max(np.array(spectra1),axis=1)[:,np.newaxis]
    max2 = np.max(np.array(spectra2),axis=1)[:,np.newaxis]

    # normed0 = np.zeros_like(max0)
    normed1 = np.zeros_like(max1)
    normed2 = np.zeros_like(max2)
    # np.true_divide(1,max0,normed0,where=(max0>0))
    np.true_divide(1,max1,normed1,where=(max1>0))
    np.true_divide(1,max2,normed2,where=(max2>0))
    # normed0 = normed0 * spectra0
    normed1 = normed1 * spectra1
    normed2 = normed2 * spectra2

    # average0 = np.mean(normed0,axis=0)
    average1 = np.mean(normed1,axis=0)
    average2 = np.mean(normed2,axis=0)

    return average1,average2




subfolders = var.list_folders


current_folder=var.current_folder



list_color=['red']

# give the bechmark spectra 

dict_files = {
    "red": f"{var.path_xml}/red_color_mito_lipid_choosenJuly22_1.xml"

}
dict_organelles = {
    "red": ["mitochondria","LD"]

}


#read the z stack raw images from different folder without unmixed
list_meta = []
for subfolder in subfolders:
    for color in list_color:
        for path_file in (Path(f"{var.path_in}")/subfolder).glob(f"zStack-{color}*.nd2"):
            print(path_file)
            dict_meta = get_nd2_size(str(path_file))
            dict_meta["folder"] = subfolder
            dict_meta["file"] = path_file.stem
            list_meta.append(dict_meta)
df_meta = pd.concat([pd.DataFrame(meta,index=[m]) for m,meta in enumerate(list_meta)])

df_meta.to_csv(f"{var.path_out}/{current_folder}/spectra/red_meta_normalized.csv",index=False)




list_benchmark = []
for color in list_color:
    df_read = read_spectra_xml(str(dict_files[color]))
    for i in range(2):
        list_benchmark.append(
            pd.DataFrame({
                "wavelength": (df_read.loc[(df_read['eType'].eq(2) & df_read['unmix'].eq(0)),'dWavelength'].to_numpy()
                              +df_read.loc[(df_read['eType'].eq(3) & df_read['unmix'].eq(0)),'dWavelength'].to_numpy()
                              )/2,
                "intensity": df_read.loc[(df_read['unmix'].eq(i) & df_read['eType'].eq(2)),"dTValue"],
                "organelle": dict_organelles[color][i]
            })
        )
df_benchmark = pd.concat(list_benchmark)
df_benchmark_max = df_benchmark.groupby("organelle").max()
df_benchmark.set_index(["organelle",'wavelength'],inplace=True)
df_benchmark["norm"] = df_benchmark.loc[:,"intensity"]/df_benchmark_max.loc[:,'intensity']
df_benchmark.reset_index(inplace=True)
##
df_meta = pd.read_csv(f"{var.path_out}/{current_folder}/spectra/red_meta_normalized.csv")
list_df = []
for folder,stem in zip(df_meta['folder'],df_meta['file']):
    # color = 'blue' 
    # print(color)
    
    file_raw = Path(f"{var.path_in}{folder}/{stem}.nd2")
    print(file_raw)
    
    file_label1 = f"{var.path_out}/{current_folder}/{var.output_folders['ilastk_segmentated']}/{dict_organelles[color][0]}_{file_raw.stem.partition('-')[2]}.tiff"
    
    print(file_label1)
    
    
    file_label2 = f"{var.path_out}/{current_folder}/{var.output_folders['ilastk_segmentated']}/{dict_organelles[color][1]}_{file_raw.stem.partition('-')[2]}.tiff"
    print(file_label2)
    
    try:
        mean1,mean2 = get_spectra_img(
                                  str(file_label1),
                                   str(file_label2),
                                  str(file_raw)
                                 )
    except FileNotFoundError:
        print(f"File Not Found: {folder}/{stem}")
        continue
    list_df.append(
        pd.DataFrame({
            "organelle":     (orga:=dict_organelles[color][0]),
            "wavelength":    df_benchmark.loc[df_benchmark['organelle'].eq(orga),'wavelength'].to_list()[:len(mean1)],
            "intensity":     mean1,
            "experiment":    folder,
            "field_of_view": stem
        })        
    )
    list_df.append(
        pd.DataFrame({
            "organelle":     (orga:=dict_organelles[color][1]),
            "wavelength":    df_benchmark.loc[df_benchmark['organelle'].eq(orga),'wavelength'].to_list()[:len(mean1)],
            "intensity":     mean2,
            "experiment":    folder,
            "field_of_view": stem
        })        
    )
df_img = pd.concat(list_df)
df_img.to_csv(f"{var.path_out}/{current_folder}/spectra/red_spectra_images.csv",index=False)


# Plot spectra:

df_img = pd.read_csv(f"{var.path_out}/{current_folder}/spectra/red_spectra_images.csv")
df_exp = df_img.groupby(['experiment','organelle','wavelength']).mean()
df_exp.reset_index(inplace=True)

df_max = df_exp.groupby(['experiment','organelle']).max()
df_max.reset_index(inplace=True)

df_exp.set_index(['experiment','organelle'],inplace=True)
df_max.set_index(['experiment','organelle'],inplace=True)

df_exp['norm'] = df_exp.loc[:,'intensity']/df_max.loc[:,'intensity']
df_exp.reset_index(inplace=True)

# color = 'blue'
for color in list_color:
    fig = go.Figure()
    for experiment in df_exp['experiment'].unique():
        if experiment == "EYrainbow_leucine_large":
            continue
        fig.add_trace(
            go.Scatter(
                x = df_exp.loc[(df_exp["experiment"].eq(experiment) & df_exp["organelle"].eq(dict_organelles[color][0])),'wavelength'],
                y = df_exp.loc[(df_exp["experiment"].eq(experiment) & df_exp["organelle"].eq(dict_organelles[color][0])),'norm'],
                name = f"{dict_organelles[color][0]}: {experiment}",
                mode="lines+markers",
                line=dict(dash="dash",shape="spline",color='green')
            )
        )
        
        
                   
        fig.add_trace(
            go.Scatter(
                x = df_benchmark.loc[(df_benchmark["organelle"].eq(dict_organelles[color][0])),'wavelength'],
                y = df_benchmark.loc[(df_benchmark["organelle"].eq(dict_organelles[color][0])),'norm'],
                name = dict_organelles[color][0],
                mode="lines+markers",
                line=dict(dash="solid",shape="spline",color='blue')
            )
        )
        
    
        
        
        
        fig.add_trace(
            go.Scatter(
                x = df_exp.loc[(df_exp["experiment"].eq(experiment) & df_exp["organelle"].eq(dict_organelles[color][1])),'wavelength'],
                y = df_exp.loc[(df_exp["experiment"].eq(experiment) & df_exp["organelle"].eq(dict_organelles[color][1])),'norm'],
                name = f"{dict_organelles[color][1]}: {experiment}",
                mode="lines+markers",
                line=dict(dash="dash",shape="spline",color='grey')
            )
        )
        
 
    
    
    fig.add_trace(
        go.Scatter(
            x = df_benchmark.loc[(df_benchmark["organelle"].eq(dict_organelles[color][1])),'wavelength'],
            y = df_benchmark.loc[(df_benchmark["organelle"].eq(dict_organelles[color][1])),'norm'],
            name = dict_organelles[color][1],
            mode="lines+markers",
            line=dict(dash="solid",shape="spline",color='red')
        )
    )
    fig.update_layout(template="simple_white")
    fig.write_html(f"{var.path_out}/{current_folder}/spectra/red_white_unmix_comparison_{color}.html")