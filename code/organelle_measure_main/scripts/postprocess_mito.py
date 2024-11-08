import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,measure
from organelle_measure.tools import batch_apply
import organelle_measure.vars_allround1data as var

def postprocess_mito(path_in,path_out):
    with h5py.File(str(path_in),'r') as f_in:
        img_in = f_in["exported_data"][:]
    img_in = (img_in[1]>img_in[0])
    img_out = measure.label(img_in)
    io.imsave(
        str(path_out),
        util.img_as_uint(img_out)
    )
    return None

list_i = []
list_o = []
for folder in var.list_folders:
    for path_binary in (Path(f"{var.path_out}")/folder/f"{var.output_folders['ilastk_segmentated']}").glob("mitochond*.h5"):
        
        print(path_binary)
        
    # for path_binary in (Path(folder_i)/folder).glob("probability_mito*.h5"):
        # (Path(f"{var.path_out}")/folder)/f"{var.output_folders['label_organelle']}"/f"label-vacuole_{path_binary.stem.partition('zStack-')[2]}.tiff"
        path_output = (Path(f"{var.path_out}")/folder)/f"{var.output_folders['label_organelle']}"/f"label-{path_binary.stem.partition('zStack-')[0]}.tiff"
        list_i.append(path_binary)
        list_o.append(path_output)
args = pd.DataFrame({
    "path_in":  list_i,
    "path_out": list_o
})

batch_apply(postprocess_mito,args)