import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
from organelle_measure.tools import open_organelles,neighbor_mean,batch_apply
# from organelle_measure.vars_allround1data import list_folders
import organelle_measure.vars_allround1data as var

def preprocess_blue(path_in,path_out,organelle):
    img_raw   = open_organelles[organelle](str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.75,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

list_in   = []
list_out  = []
list_orga = []

for folder in var.list_folders:
    # print(Path(f"{var.path_in}")/f"{folder}/unmixed_{folder}")
    for path_in in (Path(f"{var.path_in}")/f"{folder}").glob("unmixed_zStack-blue*.nd2"):
        print(path_in)
        path_peroxisome = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'peroxisome_{path_in.stem.partition("-")[2]}.tif'
        list_in.append(path_in)
        list_out.append(path_peroxisome)
        list_orga.append("peroxisome")

        path_vacuole = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'vacuole_{path_in.stem.partition("-")[2]}.tif'
        list_in.append(path_in)
        list_out.append(path_vacuole)
        list_orga.append("vacuole")
args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_blue,args)
