import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
from organelle_measure.tools import open_organelles,neighbor_mean,batch_apply
# from organelle_measure.vars_allround1data import list_folders
import organelle_measure.vars_allround1data as var

def preprocess_green(path_in,path_out,organelle):
    img_raw   = open_organelles[organelle](str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.3,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

list_in   = []
list_out  = []
list_orga = []
for folder in var.list_folders:
   
    for path_in in (Path(f"{var.path_in}")/f"{folder}").glob("unmixed_zStack-green*.nd2"):
        print(path_in)
        path_ER = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'ER_{path_in.stem.partition("-")[2]}.tif'
        print(path_ER)
        list_in.append(path_in)
        list_out.append(path_ER)
        list_orga.append("ER")

args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_green,args)


# BEGIN of ND2reader wrapper:
from nd2reader import ND2Reader

def get_nd2_size(path:str):
    with ND2Reader(path) as images:
        size = images.sizes
    return size

def load_nd2_plane(path:str,frame:str='cyx',axes:str='tz',idx:int=0):
    """read an image flexibly with ND2Reader."""
    with ND2Reader(path) as images:
        images.bundle_axes = frame
        images.iter_axes = axes
        img = images[idx]
    return img.squeeze()

# load_nd2_plane(str(path_in),frame='cyx',axes='tz',idx=0)
# get_nd2_size(str(Path(f"{var.path_in}")/f"{folder}"/"unmixed_zStack-green_whi5_ndz_33_time_25_field_1.nd2"))













