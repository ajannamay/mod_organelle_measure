import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,measure
from organelle_measure.tools import batch_apply
import organelle_measure.vars_allround1data as var
from time import sleep
from tqdm import tqdm
from termcolor import colored

def parse_meta_organelle(name):
    """name is the stem of the ORGANELLE label image file."""
    # name='label-peroxisome_blue_whi5_chx_0.1_rapa_2_time_3h_field_1.tiff'
    
    # name='ER_green_whi5_hu_2_time_75_field_1.tiff'
    # name="vacuole_blue_whi5_ndz_33_time_0_field_1.csv"
    
    organelle = name.partition("-")[2].partition("_")[0]
    if "bfa" in name:
        experiment = "bfa"
    if "hu" in name:
        experiment = "hu"
    elif "cycloheximide" in name:
        experiment = "cycloheximide"
    else:
        experiment = name.partition("whi5_")[2].partition("_")[0]
        # name.partition("cycloheximide_")[2].partition("_")[0]
    condition = name.partition(f"{experiment}_")[2].partition("_")[0]
    time = name.partition('time_')[2].partition("_")[0]
    
    field = name.partition("field_")[2]    
    return {
        "experiment": experiment,
        "condition":  condition,
        "hour":       time,
        "field":      field,
        "organelle":  organelle
    }

def measure1organelle(path_in,path_cell,path_out):
    # parse metadata from filename
    name = Path(path_in).stem
    meta = parse_meta_organelle(name)

    img_orga = io.imread(str(path_in))
    
    img_cell = io.imread(str(path_cell))
    
    dfs = []
    for cell in measure.regionprops(img_cell):
        meta["idx-cell"] = cell.label
        min_row, min_col, max_row, max_col = cell.bbox
        img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col]
        img_cell_crop = cell.image
        for z in range(img_orga_crop.shape[0]):
            img_orga_crop[z] = img_orga_crop[z]*img_cell_crop
        if not meta["organelle"] == "vacuole":
            
            measured_orga = measure.regionprops_table(
                img_orga_crop,
                properties=('label','area','bbox_area','bbox',)
            )
        else:
            vacuole_area = 0
            vacuole_bbox_area = 0
            bbox0,bbox1,bbox2,bbox3,bbox4,bbox5 = 0,0,0,0,0,0
            for z in range(img_orga_crop.shape[0]):
                vacuole = measure.regionprops_table(
                    img_orga_crop[z],
                    properties=('label','area','bbox_area','bbox',)
                )
                if len(vacuole["area"]) == 0:
                    continue
                if (maxblob:=max(vacuole["area"])) > vacuole_area:
                    vacuole_area = maxblob
                    idxblob = np.argmax(vacuole["area"])
                    vacuole_bbox_area = vacuole["bbox_area"][idxblob]
                    bbox0,bbox3 = z,z
                    bbox1,bbox2,bbox4,bbox5 = [vacuole[f"bbox-{i}"][idxblob] for i in range(4)]
            if vacuole_area==0:
                continue
            measured_orga = {
                'label': [0],
                'area':  [vacuole_area],
                "bbox_area": [vacuole_bbox_area],
                "bbox-0": [bbox0],
                "bbox-1": [bbox1],
                "bbox-2": [bbox2],
                "bbox-3": [bbox3],
                "bbox-4": [bbox4],
                "bbox-5": [bbox5],
                # 'intensity_mean':[0]
               
            }
        result = meta | measured_orga
        dfs.append(pd.DataFrame(result))
    if len(dfs) == 0:
        print(f">>> {path_out} has no cells, skipped.")
        return None
    df_orga = pd.concat(dfs,ignore_index=True)
    df_orga.rename(columns={'label':'idx-orga',"area":"volume-pixel",'bbox_area':'volume-bbox','intensity_mean':'intensity_mean'},inplace=True)
    df_orga.to_csv(str(path_out),index=False)
    print(colored(f">>> finished {path_out.stem}.","red",attrs=["bold"]))
    return None



organelles= {
    "peroxisome":"blue",
    "vacuole":"blue",
    "ER":"green",
    "golgi":"yellow",
    "mitochondria":"red",
    "LD":"red"
    
}


list_i = []
list_c = []
list_o = []
err_i  = []


for folder in var.list_folders:
    for path_cell in (Path(f"{var.path_out}")/folder/f"{var.output_folders['segmented_cell']}").glob("*.tif"):
        # print(path_cell)
        
        for org in tqdm(organelles, desc='Processing organelles'):
           
            print("====",colored("Organelle:","magenta",attrs=["bold"])
                  ,colored(f"{org}","red",attrs=["bold"]),"===========")
            
            sleep(0.5)
            
            
            path_i = Path(f"{var.path_out}")/folder/f"{var.output_folders['label_organelle']}"/f"label-{org}_{organelles[org]}_whi5_{path_cell.stem.partition('_')[2]}.tiff"
            
            path_o = Path(f"{var.path_out}")/folder/f"{var.output_folders['measurements']}"/f"{org}_{path_i.stem.partition('_')[2]}.csv"
            
            # print(path_o)
            if not path_i.exists():
                err_i.append(path_i)

            list_i.append(path_i)
            list_c.append(path_cell)
            list_o.append(path_o)
        
args = pd.DataFrame({
    "path_in":   list_i,
    "path_cell": list_c,
    "path_out":  list_o
})

batch_apply(measure1organelle,args)



