import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util
from organelle_measure.tools import skeletonize_zbyz,batch_apply
import organelle_measure.vars_allround1data as var

def postprocess_ER(path_in,path_out):
    with h5py.File(str(path_in)) as f_in:
        img_in = f_in["exported_data"][:]
    img_in = (img_in[1]>img_in[0])
    img_ske = skeletonize_zbyz(img_in)
    io.imsave(
        str(path_out),
        util.img_as_ubyte(img_ske)
    )
    return None

folders = [
    "EYrainbow_glucose",
    "EYrainbow_glucose_largerBF",
    "EYrainbow_rapamycin_1stTry",
    "EYrainbow_rapamycin_CheckBistability",
    "EYrainbow_1nmpp1_1st",
    "EYrainbowWhi5Up_betaEstrodiol",
    "EYrainbow_leucine_large",
    "EYrainbow_leucine"
]
folder_i = "./images/preprocessed/"
folder_o = "./images/labelled/"

list_i = []
list_o = []
for folder in var.list_folders:
    for path_binary in (Path(f"{var.path_out}")/folder/f"{var.output_folders['ilastk_segmentated']}").glob("ER*.h5"):
        
        print(path_binary)
        path_output = (Path(f"{var.path_out}")/folder/f"{var.output_folders['label_organelle']}")/f"label-ER_{path_binary.stem.partition('_')[2]}.tiff"
        list_i.append(path_binary)
        list_o.append(path_output)
args = pd.DataFrame({
    "path_in":  list_i,
    "path_out": list_o
})

batch_apply(postprocess_ER,args)
# 