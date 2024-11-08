
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import segmentation,measure,io,util
from organelle_measure.yeaz import yeaz_preprocesses,yeaz_label
from organelle_measure.tools import load_nd2_plane,batch_apply
import organelle_measure.vars_allround1data as var
# import scripts.folder_details as fd
# from organelle_measure.vars_allround1data import list_folders

def segment_cells(path_in,path_out):
    img_i = load_nd2_plane(str(path_in),frame='yx',axes='t',idx=0)
    for prep in yeaz_preprocesses:
        img_i = prep(img_i)
    img_b = yeaz_label(img_i,min_dist=4)
    img_b = segmentation.clear_border(img_b)
    properties = measure.regionprops(img_b)
    for prop in properties:
        if prop.area < 50: # hard coded threshold, bad
            img_b[img_b==prop.label] = 0
    img_b = measure.label(img_b)
    img_o = np.zeros((512,512),dtype=int) # hard coded size, bad
    shape0,shape1 = img_b.shape
    img_o[:shape0,:shape1] = img_b

    io.imsave(str(path_out),util.img_as_uint(img_o))
    print(f"...{path_out}")
    return None

list_in = []
list_out = []

for folder in var.list_folders:
    for file_cell in Path(f"{var.path_in}/{folder}").glob("camera*.nd2"):
        print(file_cell)
        list_in.append(file_cell)
        file_segm = Path(f"{var.path_out}/{folder}/{var.output_folders['segmented_cell']}")/f"CellSeg-{file_cell.stem.partition('_')[2]}.tif"
        list_out.append(file_segm)
        print("output_pat",file_segm)
        
args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out
})


# print(list_in)
batch_apply(segment_cells,args)

