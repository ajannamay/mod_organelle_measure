import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,segmentation
from organelle_measure.tools import batch_apply
import organelle_measure.vars_allround1data as var
from termcolor import colored
from time import sleep
from tqdm import tqdm
import os
os.environ['Tf_CPP_MIN_LOG_LEVEL']='2'


def postprocess_globular(path_in,path_ref,path_out):
    with h5py.File(str(path_in),'r') as f_in:
        img_in = f_in["exported_data"][:]
    img_in = (img_in[1]>img_in[0])
    img_ref = io.imread(str(path_ref))
    
    # idx_max = feature.peak_local_max(img_ref,min_distance=1)
    # img_max = np.zeros_like(img_ref,dtype=bool)
    # img_max[tuple(idx_max.T)] = True

    img_out = segmentation.watershed(-img_ref,mask=img_in)
    io.imsave(
        str(path_out),
        util.img_as_uint(img_out)
    )
    return None


organelles = [
    "peroxisome",   
    "LD",
    "golgi"
]



list_i   = []
list_ref = []
list_o   = []

for folder in var.list_folders:
    for organelle in tqdm(organelles, desc='Processing organelles'):
        # tqdm(range(len(organelles)))
        sleep(2)
        
        for path_binary in (Path(f"{var.path_out}")/folder/f"{var.output_folders['ilastk_segmentated']}").glob(f"{organelle}*.h5"):
            
            # print(path_binary)

            path_output = (Path(f"{var.path_out}")/folder)/f"{var.output_folders['label_organelle']}"/f"label-{path_binary.stem.partition('-')[0]}.tiff"
            
            path_ref = (Path(f"{var.path_out}")/folder)/f"{var.output_folders['ilastk_segmentated']}"/f"{path_binary.stem.partition('-_')[0]}.tiff"
            # print(path_ref)
            list_i.append(path_binary)
            list_ref.append(path_ref)
            list_o.append(path_output)
args = pd.DataFrame({
    "path_in":  list_i,
    "path_ref": list_ref,
    "path_out": list_o
})

batch_apply(postprocess_globular,args)