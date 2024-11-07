import pandas as pd
from pathlib import Path
from skimage import io,measure
from organelle_measure.tools import batch_apply
# from organelle_measure.vars_allround1data import list_folders
import organelle_measure.vars_allround1data as var

def parse_meta_cell(name):
    """name is the stem of the ORGANELLE label image file."""
    
    # name="vacuole_blue_whi5_ndz_33_time_0_field_1.csv"
    if "bfa" in name:
        experiment = "bfa"
    if "hu" in name:
        experiment = "hu"
    elif "cycloheximide" in name:
        experiment = "cycloheximide"
    else:
        experiment =  name.partition("whi5_")[2].partition("_")[0]
    condition = name.partition(f"{experiment}_")[2].partition("_")[0]
    time = name.partition('time_')[2].partition("_")[0]
    field = name.partition("field_")[2]    
    return {
        "experiment": experiment,
        "condition":  condition,
        "hour":       time,
        "field":      field,
    }

def measure1cell(path_in,path_out):
    img_cell = io.imread(str(path_in))
    name = Path(path_in).stem
    meta = parse_meta_cell(name)
    measured = measure.regionprops_table(
                    img_cell,
                    properties=('label','area','centroid','bbox','eccentricity')
               )
    result = meta | measured
    # meta.update(measured)
    # result = meta
    
    df = pd.DataFrame(result)
    df.rename(columns={'label':'idx-cell'},inplace=True)
    df.to_csv(str(path_out),index=False)
    return None


list_i = []
list_o = []




for folder in var.list_folders:
    for path_i in (Path(f"{var.path_out}")/folder/f"{var.output_folders['segmented_cell']}").glob("*.tif"):
        # print(path_i)
        path_o = (Path(f"{var.path_out}")/folder/f"{var.output_folders['measurements']}")/f"cell_{path_i.stem.partition('-')[2]}.csv"
        
        print(path_o)
        list_i.append(path_i)
        list_o.append(path_o)

args = pd.DataFrame({
    "path_in":   list_i,
    "path_out":  list_o
})

batch_apply(measure1cell,args)


