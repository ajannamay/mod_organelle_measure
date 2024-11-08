

from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog
from os.path import expanduser

def choose_directory():
        print("Hello1")
        input_dir = QFileDialog.getExistingDirectory(None, 'Select a folder:', expanduser("~"))
     
        return input_dir



list_folders = [
    "Aug07"
]


current_folder=list_folders[0]



####new strain

# path_measure="../../../Experiments/2024/new_strain"
# path_in="../../../Experiments/2024/new_strain/"
# path_out="../../../Analysis/2024/new_strain"


# exec(open('C:Users/admin/GitHub/mod_organelle_measure/code/organelle_measure_main/scripts/_runFirst.py').read())




# # ####noise_tradeoff

# path_measure="C:/Users/admin/GitHub/mod_organelle_measure/experiments"
path_in="../../experiments"
path_out="../../Analysis"


output_folders={
    "segmented_cell":"segmented_cell",
    "organelle_unmixed":"organelle_unmixed",
    "ilastk_segmentated":"ilastk_segmentated",
    "label_organelle":"label_organelle",
    "spectra":"spectra",
    "measurements":"measurements",
    "ilastk_control_libfile":"ilastk_control_libfile"

    }



####CU

# path_measure="../../../Experiments/2024/CU"
# path_in="../../../Experiments/2024/CU/"
# path_out="../../../Analysis/2024/CU"



####chxrapa

# path_measure="../../../Experiments/2024/chxrapa"
# path_in="../../../Experiments/2024/chxrapa/"
# path_out="../../../Analysis/2024/chxrapa"






# out_fig1="../../../Analysis/2023_analysis/csv2_figure/nov"

# #for reading files
# out_chxrapa="../../../Analysis/2024/chxrapa"
# out_hu="../../../Analysis/2023_analysis/HU"
# out_chx="../../../Analysis/2023_analysis/CHX"
# out_bfa="../../../Analysis/2023_analysis/BFA"
# out_control="../../../Analysis/2023_analysis/control"
# chx_im="../../../Analysis/2023_analysis/CHX/chx_images"


# ###==========xml file ie benchmark spectra
# path_xml="../../../Experiments/2023/xml_fileJan31"

####chx

# path_measure="../../../Experiments/2023/CHX"
# path_in="../../../Experiments/2023/CHX/"
# path_out="../../../Analysis/2023_analysis/CHX"

### BFA
# path_measure="../../../Experiments/2023/BFA"
# path_in="../../../Experiments/2023/BFA/"
# path_out="../../../Analysis/2023_analysis/BFA"


####HU

# path_measure="../../../Experiments/2023/HU"
# path_in="../../../Experiments/2023/HU/"
# path_out="../../../Analysis/2023_analysis/HU"






##measurements

# subfolders_chxrapa = [
    
#    "April03",
  
#   "April04",
#   "April16",
  
#     ]


# # path_measure="../../../Analysis/2023_analysis/CHX/chx_measurements"



# growth_rate_chx={

#     "0":"0.5",
#     "0.01":"0.4",
#     "0.1":"0.29",
#     "0.5":"0.16",
#     "2":"0.13",

#     }




# growth_rate_hu={

#     "0":"0.5",
#     "1.2":"0.44",
#     "4":"0.36",
 

#     }



# raw2sp=[
#     # "July22",
#     # "July25",
#     # "July28",
#     # "Aug01",
#     # "Aug02",
#     # "Aug05",
#     # "Aug06",
#     # "Aug07",
# # "Aug17",
# "Aug19",
# "Aug20",
# "Aug21",
# "Sept07",
# "sept29",
# "Sept30",
# "Oct02",
# "Oct03",
# # "Oct04",
# # "Oct05",
# "Oct07",
# ]



# subfolders_chx = [
    
#     "Aug01",
#     "Aug02",
#     "Aug07",
    
#     #========================0.01================#
#     "Aug05",
#     "Aug19",
#     "Sept30",
#     "Oct07",
#     #========================0.1================#
#     "Aug06",
#     "Aug20",
#     "Oct02",
#     # "Oct04",
#     #========================0.5================#
#     # "Aug17",
#     "Aug21",
#     # "Oct05",
#     #========================2.0================#
#     "Sept07",
#     "Sept29",
#     "Oct03",
    
# ]
# subfolders_bfa = [
#     #========================25================#
    
#     "Nov27","Nov28","Sept11",
#     #========================75================#
#     "Nov30","Sept12",
#     #========================100================#
    
#     "Sept13",
    
# ]

# subfolders_control = [
   
#     "Aug01",
#     "Aug02",
#     "Aug07",
# ]


# subfolders_hu = [
#     #========================1.2================#
#     "Dec12",
#     "Nov06",
#     "Nov08",
#     #========================4================#
#     "Oct21"
    
    
# ]


# organelles={
    
    
#     "peroxisome":"blue",
#     "vacuole":"blue",
#     "ER":"green",
#     "golgi":"yellow",
#     "mitochondria":"red",
#     "LD":"red",
#     "non-organelle":"non-organelle"
    
    
#     }






# experiments = {
#     "control":"control",
#     "cycloheximide": "Cycloheximide",
#     "bfa":     "bfa",
#     "hu":"hu"
# }

# growth_rate={
    
#     "growth_chx":{"0.0":0.5,
#         "0.01":0.4,
#                   "0.1":0.29,
#                   "0.5":0.16,
#                   "2.0":0.13},
    
    
#     # "growth_chx":[0.4,0.29,0.16,0.13],
#     "growth_bfa":{
    
#         "25.0":0.4,
#         "75.0":0.29,
#         "100.0":0.13
        
#         },
#     "growth_hu":
        
#         {
#             "1.2":0.4,
                
#                 "4":0.29
            
            
#             }
        
    
    
#     }




















# from statsmodels.stats.weightstats import ttest_ind
# from statannotations.Annotator import Annotator
# from scipy.stats import bootstrap
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
# from matplotlib import container
# import scipy
# from sklearn import metrics

# ##get the number of cell in each condition
# for cond in df_condition:
#     df_bycell_volume=df_bycell.loc[ df_bycell["condition"].eq(cond)]
#     print("condition:",cond,"independent measurement:",df_bycell_volume["folder"].unique())
#     # print(df_bycell_volume1)
#     print("Number of cells in each condition",cond,df_bycell_volume.shape[0]/7)

'''

#read the files
df_bycell_control,df_orga_all = read_results(Path(out_control),subfolders_control,(px_x,px_y,px_z))
df_bycell_chx,df_orga_all = read_results(Path(out_chx),subfolders_chx,(px_x,px_y,px_z))
df_bycell_bfa,df_orga_all = read_results(Path(out_bfa),subfolders_bfa,(px_x,px_y,px_z))
df_bycell=pd.concat([df_bycell_control,df_bycell_chx,df_bycell_bfa])

# Apply the custom function to create the new column
df_bycell['experiment'] = df_bycell.apply(custom_function, axis=1)

'''


    # cmap = sns.diverging_palette(220, 10, as_cmap=True)
    #  # Draw the heatmap with the mask and correct aspect ratio
    # # More details at https://seaborn.pydata.org/generated/seaborn.heatmap.html
    # sns.heatmap(
    #     corrcoef,          # The data to plot
    #     # mask=mask,     # Mask some cells
    #     cmap=cmap,     # What colors to plot the heatmap as
    #     annot=True,    # Should the values be plotted in the cells?
    #     vmax=.3,       # The maximum value of the legend. All higher vals will be same color
    #     vmin=-.3,      # The minimum value of the legend. All lower vals will be same color
    #     center=0,      # The center value of the legend. With divergent cmap, where white is
    #     square=True,   # Force cells to be square
    #     linewidths=.5, # Width of lines that divide cells
    #     cbar_kws={"shrink": .5}  # Extra kwargs for the legend; in this case, shrink by 50%
    # )
    
    
    

# from termcolor import colored

# # Sample variables
# org = "Golgi"
# var = "mean"
# value = 42

# # Define colors
# org_color = 'red'
# var_color = 'blue'
# value_color = 'green'

# # Print colored text with variable values
# print(colored(f"Organelle: {org}", org_color))
# print(colored(f"Variable: {var}", var_color))
# print(colored(f"Value: {value}", value_color))