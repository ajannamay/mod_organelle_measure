import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,morphology,measure
from organelle_measure.tools import skeletonize_zbyz,watershed_zbyz
from organelle_measure.tools import find_complete_rings,better_vacuole_img,batch_apply
import organelle_measure.vars_allround1data as var

# path_in=args["path_in"][0]
# path_out=args["path_out"][0]
# path_cell=args["path_cell"][0]

def postprocess_vacuole(path_in,path_cell,path_out):
    with h5py.File(str(path_in),'r') as f_in:
        img_orga = f_in["exported_data"][:]
    img_orga = np.argmax(img_orga,axis=0)
    img_orga = (img_orga>0)

    img_cell = io.imread(str(path_cell))

    img_skeleton  = skeletonize_zbyz(img_orga)

    img_core      = find_complete_rings(img_skeleton)
    
    # img_vacuole   = better_vacuole_img(img_core,img_watershed)
    img_vacuole = np.zeros_like(img_core,dtype=int)
    for z in range(img_vacuole.shape[0]):
        sample = img_core[z]
        candidates = np.unique(sample[img_cell>0])
        for color in candidates:
            if len(np.unique(img_cell[sample==color]))==1:
                img_vacuole[z,sample==color] = color

    io.imsave(
        str(path_out),
        util.img_as_uint(img_vacuole) 
    )
    return None



list_i = []
list_c = []
list_o = []
for folder in var.list_folders:
    for path_cell in (Path(f"{var.path_out}")/folder/
                      f"{var.output_folders['segmented_cell']}").glob("*.tif"):
        print(path_cell)
        path_binary = (Path(f"{var.path_out}")/folder/
                       f"{var.output_folders['ilastk_segmentated']}")/f"vacuole_blue_{path_cell.stem.partition('-')[2]}.h5"
        print(path_binary)
        path_output = (Path(f"{var.path_out}")/folder)/f"{var.output_folders['label_organelle']}"/f"label-{path_binary.stem.partition('-_')[0]}.tiff"
        
        list_i.append(path_binary)
        list_c.append(path_cell)
        list_o.append(path_output)

args = pd.DataFrame({
    "path_in":   list_i,
    "path_cell": list_c,
    "path_out":  list_o
})

batch_apply(postprocess_vacuole,args)
