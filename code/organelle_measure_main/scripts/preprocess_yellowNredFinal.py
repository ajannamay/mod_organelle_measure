import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
from organelle_measure.tools import open_organelles,neighbor_mean,batch_apply
# from organelle_measure.vars_allround1data import list_folders
import organelle_measure.vars_allround1data as var

# clean the yellow and red channels
def preprocess_yellowNred(path_in,path_cell,path_out,organelle):
    img_cell  = io.imread(str(path_cell))
    img_raw   = open_organelles[organelle](str(path_in))
    # # special for leucine large unmixed-red
    # img = io.imread(str(path_in))
    # img_raw = np.zeros([int(img.shape[0]/2),*img.shape[1:]],dtype=int)
    # if organelle == "mitochondria":
    #     for z in range(img_raw.shape[0]):
    #         img_raw[z] = img[2*z]
    # elif organelle == "LD":
    #     for z in range(img_raw.shape[0]):
    #         img_raw[z] = img[2*z+1]
    img_bkgd  = neighbor_mean(img_raw,img_cell)
    img_clean = img_raw - img_bkgd
    img_clean[img_clean<0] = 0
    img_gaussian = filters.gaussian(img_clean,sigma=0.75,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

list_in   = []
list_cell = []
list_out  = []
list_orga = []
for folder in var.list_folders:
    for path_cell in (Path(f"{var.path_out}")/f"{folder}"/f"{var.output_folders['segmented_cell']}").glob("Cell*.tif"):
        print(path_cell)
        
        # (Path(f"{var.path_in}")/f"{folder}/unmixed_{folder}")
        # unmixed_zStack-yellow_whi5_control_0_time_0_field_1.nd2
        
        path_yellow = Path(f"{var.path_in}")/f"{folder}"/f"unmixed_zStack-yellow_{path_cell.stem.partition('-')[2]}.nd2"
        print(path_yellow)
        
        
        path_yellow = Path(f"{var.path_in}")/f"{folder}"/f"unmixed_zStack-yellow_{path_cell.stem.partition('-')[2]}.nd2"
        print(path_yellow)
        
        path_golgi = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'golgi_{path_yellow.stem.partition("zStack-")[2]}.tif'
        print(path_golgi)
        list_in.append(path_yellow)
        list_cell.append(path_cell)
        list_out.append(path_golgi)
        list_orga.append("golgi")

        if folder == "EYrainbow_leucine":
            path_red = Path("./data/raw")/folder/f"unmixed-red_{path_cell.stem.partition('_')[2]}.ome.tif"
        else:
            path_red  = Path(f"{var.path_in}")/f"{folder}"/f"unmixed_zStack-red_{path_cell.stem.partition('-')[2]}.nd2"
        path_mitochondria = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'mitochondria_{path_red.stem.partition("zStack-")[2]}.tif'

        list_in.append(path_red)
        list_cell.append(path_cell)
        list_out.append(path_mitochondria)
        list_orga.append("mitochondria")

        path_LD = Path(f"{var.path_out}")/folder/f"{var.output_folders['organelle_unmixed']}"/f'LD_{path_red.stem.partition("zStack-")[2]}.tif'

        list_in.append(path_red)
        list_cell.append(path_cell)
        list_out.append(path_LD)
        list_orga.append("LD")   
args = pd.DataFrame({
    "path_in": list_in,
    "path_cell": list_cell,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_yellowNred,args)

