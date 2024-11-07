# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:12:47 2023

@author: Anang
"""


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from organelle_measure.data import read_results
import organelle_measure.vars_allround1data as var
import scipy.stats as stats
from termcolor import colored
import organelle_measure.an_functions as an
from organelle_measure.data import read_results_pixel



###data frame adding cell_all=pd.concat([cell_bfa,cell_chx,cell_control])
font_size=12
chx_im=f"{var.chx_im}"
px_x,px_y,px_z = 0.41,0.41,0.20

##
out_hu=f"{var.out_hu}"
out_chx=f"{var.out_chx}"
out_bfa=f"{var.out_bfa}"
out_control=f"{var.out_control}"
out_fig="../../../Analysis/2023_analysis/csv2_figure/2023"
#####get the folder and variables 
subfolders_chx=var.subfolders_chx
subfolders_bfa=var.subfolders_bfa
subfolders_control=var.subfolders_control
organelles=var.organelles
###=====================================================
subfolders_hu = [
    "Oct09",
    "Oct11",
    "Oct20",

]

experiments = {
    "control":"control",
    "cycloheximide": "Cycloheximide",
    "bfa":     "bfa",
    "hu":"hu"
}

exp_names = experiments.keys()
exp_names = list(exp_names)
exp_folder = [experiments[i] for i in exp_names]
#extremes for each cases
extremes = {
    "control":[0.],
     "cycloheximide": [ 0.01,0.1,0.5,2],
     "bfa":           [ 25,75, 100.],
     "hu":[0.2,2,4]
     
}

df_condition=extremes["control"]+extremes["cycloheximide"]+extremes["bfa"]+extremes["hu"]
folders=[subfolders_hu,subfolders_bfa,subfolders_chx,subfolders_control]

list_colors = {
    "cycloheximide":     [1,2,3,4,0,5],
    "bfa":     [1,2,3,4,0],
    "control":   [1,0],

}


for exp in extremes:
    for values in extremes[f"{exp}"]:
        print("experiment",exp,"values",values)
        
        
##hu
df_bycell_hu = read_results(Path(out_hu),subfolders_hu,(px_x,px_y,px_z))

df_bycell_hu1,df_orga_all1 = read_results_pixel(Path(out_hu),subfolders_hu)



df_bycell_chx=pd.read_csv("../../../Analysis/Analysis_files/df_bycell_Oct10_2023.csv")


xx=df_bycell_chx[df_bycell_chx["folder"] != "Sept29"]

subfolders=["Sept29"]
df_bycell_ch05,df_orga_all = read_results(Path(out_chx),subfolders,(px_x,px_y,px_z))

df_bycell_ch05['experiment'] = df_bycell_ch05.apply(an.custom_function, axis=1)

df_bycell_concated=pd.concat([df_bycell_ch05,df_bycell_chx])


var="mean"
exp_names=["cycloheximide"]
an.get_each_date_field(organelles,exp_names,df_bycell_ch05,var)




df_bycell_concated=pd.concat([df_bycell_hu,df_bycell_chx])

df_bycell2=df_bycell_concated[df_bycell_concated["cell-volume"]<1000]

df_bycell2=df_bycell2[df_bycell2["cell-volume"]>50]
 #eleminating very high volume cells
df_bycell=df_bycell2

df_bycell=df_bycell_concated

df_bycell["cell-volume"].min()
# DATAFRAME FOR CORRELATION COEFFICIENT
pv_bycell = df_bycell.set_index(['folder','condition','field','idx-cell','experiment'])
df_corrcoef = pd.DataFrame(index=pv_bycell.loc[pv_bycell["organelle"].eq("ER")].index)
df_corrcoef.loc[:,'effective-length'] = np.sqrt(pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-area']/np.pi)
df_corrcoef.loc[:,'cell-area'] = pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-area']
df_corrcoef.loc[:,'cell-volume'] = pv_bycell.loc[pv_bycell["organelle"].eq("ER"),'cell-volume']

properties = []
for orga in [*organelles,"non-organelle"]:
    for prop in ['mean','count','total','total-fraction']:
        if (orga in ["ER","vacuole","non-organelle"]) and (prop in ["count","mean"]):
            continue
        prop_new = f"{prop}-{orga}"
        properties.append(prop_new)
        df_corrcoef.loc[:,prop_new] = pv_bycell.loc[pv_bycell["organelle"]==orga,prop]
df_corrcoef.reset_index(inplace=True)
columns = ['effective-length','cell-area','cell-volume',*properties]


#########get the number of cells
#########===========================================================###
###########=======================ploting=======================#######
#########===========================================================###
##==========================================##
###==========================================#

var="mean"

exp_names=["control"]
an.get_each_date_field(organelles,exp_names,df_bycell,var)




experiment="cycloheximide"
property="total-fraction"
an.make_pca_plots1(df_bycell,organelles,experiment,out_fig,list_colors,property,has_volume=False,is_normalized=True,non_organelle=False)

colors = np.random.randint(100, size=(100))

an.get_details(exp_names,df_bycell)
an.cell_volume(extremes,df_bycell)
an.cell_volume("cycloheximide",df_bycell)
an.corelation_plot(extremes,df_corrcoef,columns,properties,cmap=False)




exp_names=["control","cycloheximide","bfa"]
control_org=an.get_independent_measurement(organelles,exp_names,df_bycell,"total-fraction",with_error_bar=True)

group_control=control_org.loc[control_org["condition"]==0.0,"mean"]

exp_names=["cycloheximide"]
chx_org=an.get_independent_measurement(organelles[3],exp_names,df_bycell,"total-fraction",with_error_bar=True)

#################PCA
experiment="cycloheximide"
property="total-fraction"
data1=df_bycell[df_bycell["experiment"]=="control"]
data2=df_bycell[df_bycell["experiment"]==f"{experiment}"]
groups=[data1["condition"].max(),data2["condition"].max()]

folder1 =experiment# experiments[experiment]
name = f"{'all-conditions' if groups is None else 'extremes'}_{'has-cytoplasm' if non_organelle else 'no-cytoplasm'}_{'cell-volume' if has_volume else 'organelle-only'}_{property}_{'norm-mean-std' if is_normalized else 'raw'}"
print("PCA Anaysis: ",folder1,name)    
df_orga_perfolder=pd.concat([data1,data2])

df_orga_perfolder.set_index(["experiment","condition","field","idx-cell","folder"],inplace=True)
idx = df_orga_perfolder.groupby(["experiment","condition","field","idx-cell","folder"]).count().index

columns = [*organelles,"non-organelle"] if non_organelle else organelles
df_pca = pd.DataFrame(index=idx,columns=columns)
num_pc = 7 if non_organelle else 6

for orga in columns:
    print(orga)
    df_pca[orga]=df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq(orga),property]


if has_volume:
    df_pca["cell-volume"] = df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq("ER"),"cell-volume"]
    columns = ["cell-volume",*columns]
    num_pc += 1

if is_normalized:
    for col in columns:
        df_pca[col] = (df_pca[col]-df_pca[col].mean())/df_pca[col].std()
        
        
sns.pairplot(df_pca.drop(columns=[ 'field','idx-cell','folder','condition']), hue="experiment", height=3, diag_kind="kde")
plt.show()        


df_pca.reset_index(inplace=True)

np_pca = df_pca[columns].to_numpy()

sns.jointplot(x="LD", y="peroxisome", data=df_pca)
plt.show()

sns.FacetGrid(df_pca, hue="experiment", size=5) \
    .map(plt.scatter, "LD", "peroxisome") \
    .add_legend()
plt.show()


import seaborn as sns

mean = np.mean(np_pca, axis=0)
std = np.std(np_pca, axis=0)
X = (np_pca - mean) / std








z1=df_pca.condition.unique()
from sklearn.manifold import TSNE

tsne = TSNE(n_components=2, perplexity=30, n_iter=300, random_state=0)
X_tsne = tsne.fit_transform(X)

# Create a scatter plot of the t-SNE results
plt.figure(figsize=(8, 6))
colors =['navy', 'turquoise', 'darkorange',"blue","green"]
for i in range(4):
    plt.scatter(X_tsne[y == i, 0], X_tsne[y == i, 1], c=colors[i], label=z1[i])

plt.title('t-SNE Visualization of Iris Dataset')
plt.legend()
plt.show()


figproj = plt.figure(figsize=(15,12))
# ax = figproj.add_subplot(projection="3d")
for color, i, target_name in zip(colors, lb1,lb1):
    ax.scatter(X_tsne[y == i, 0], X_tsne[y == i, 1], color=color, alpha=0.8, lw=lw,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of IRIS dataset')



n_bins = 50

fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
fig, ax = plt.subplots(tight_layout=True)
ax.hist2d(X[:,0], X[:,1])
for i in range(5):
    axs.hist(X[:,i], bins=n_bins)
    
    
    counts, bins = np.histogram(X[:,i])
    plt.stairs(counts, bins)


sns.jointplot(x=X[:,0], data=X, size=5)
plt.show()


# Perform PCA
pca = PCA(n_components=4)  # Reduce to 2 principal components
X_r = pca.fit(X).transform(X)
y = df_pca.condition
# Percentage of variance explained by each of the selected components
explained_variance = pca.explained_variance_ratio_

plt.bar(range(4), explained_variance, align='center', alpha=0.5)
# plt.xticks(range(2), ['PC1', 'PC2'])

colors =['navy', 'turquoise', 'darkorange',"blue","green"]
lw = 2

lb1=[0,2]

figproj = plt.figure(figsize=(15,12))
ax = figproj.add_subplot(projection="3d")
for color, i, target_name in zip(colors, lb1,lb1):
    ax.scatter(X_r[y == i, 1], X_r[y == i, 2],X_r[y == i, 3], color=color, alpha=0.8, lw=lw,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of IRIS dataset')









# # PCA 
# def make_pca_plots(experiment,property,groups=None,has_volume=False,is_normalized=False,non_organelle=False):
#     folder = experiments[experiment]
#     name = f"{'all-conditions' if groups is None else 'extremes'}_{'has-cytoplasm' if non_organelle else 'no-cytoplasm'}_{'cell-volume' if has_volume else 'organelle-only'}_{property}_{'norm-mean-std' if is_normalized else 'raw'}"
#     print("PCA Anaysis: ",folder,name)

#     df_orga_perfolder = df_bycell[df_bycell["folder"].eq(folder)]
#     df_orga_perfolder.set_index(["condition","field","idx-cell"],inplace=True)
#     idx = df_orga_perfolder.groupby(["condition","field","idx-cell"]).count().index
    
#     columns = [*organelles,"non-organelle"] if non_organelle else organelles
#     df_pca = pd.DataFrame(index=idx,columns=columns)
#     num_pc = 7 if non_organelle else 6
    
#     for orga in columns:
#         df_pca[orga] = df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq(orga),property]
    
#     if has_volume:
#         df_pca["cell-volume"] = df_orga_perfolder.loc[df_orga_perfolder["organelle"].eq("ER"),"cell-volume"]
#         columns = ["cell-volume",*columns]
#         num_pc += 1
    
#     if is_normalized:
#         for col in columns:
#             df_pca[col] = (df_pca[col]-df_pca[col].mean())/df_pca[col].std()
    
#     df_pca.reset_index(inplace=True)

#     # Find the the direction of the condition change:
#     df_centroid = df_pca.groupby("condition")[columns].mean()
#     fitter_centroid = LinearRegression(fit_intercept=False)
#     np_centroid = df_centroid.to_numpy()
#     fitter_centroid.fit(np_centroid,np.ones(np_centroid.shape[0]))
    
#     vec_centroid_start = df_centroid.loc[groups[0],:].to_numpy()
#     vec_centroid_start[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_start[:-1]))/fitter_centroid.coef_[-1]

#     vec_centroid_ended = df_centroid.loc[groups[-1],:].to_numpy()
#     vec_centroid_ended[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_ended[:-1]))/fitter_centroid.coef_[-1]

#     vec_centroid = vec_centroid_ended - vec_centroid_start
#     vec_centroid = vec_centroid/np.linalg.norm(vec_centroid)
#     np.savetxt(f'{Path(f"{out_fig}")}/condition-vector_{folder}_{name}.txt',vec_centroid)

#     # Get Principal Components (PCs)
#     np_pca = df_pca[columns].to_numpy()
#     pca = PCA(n_components=num_pc)
#     pca.fit(np_pca)
#     pca_components = pca.components_
#     pca_var_ratios = pca.explained_variance_ratio_

#     # Calculate cos<condition,PCs>, and realign the PCs
#     cosine_pca = np.dot(pca_components,vec_centroid)
#     for c in range(len(cosine_pca)):
#         if cosine_pca[c] < 0:
#             pca_components[c] = -pca_components[c]
#     cosine_pca = np.abs(cosine_pca)
#     # save and plot PCs without sorting.
#     np.savetxt(f'{Path(f"{out_fig}")}/cosine_{folder}_{name}.txt',cosine_pca)
#     np.savetxt(f'{Path(f"{out_fig}")}/pca-components_{folder}_{name}.txt',pca_components)
#     np.savetxt(f'{Path(f"{out_fig}")}/pca-explained-ratio_{folder}_{name}.txt',pca_var_ratios)
#     fig_components = px.imshow(
#         pca_components,
#         x=columns, y=[f"PC{i}" for i in range(num_pc)],
#         color_continuous_scale="RdBu_r", color_continuous_midpoint=0
#     )
#     fig_components.write_html(f'{Path(f"{out_fig}")}/pca_data/pca_components_{folder}_{name}.html')

#     # Sort PCs according to the cosine
#     arg_cosine = np.argsort(cosine_pca)[::-1]
#     pca_components_sorted = pca_components[arg_cosine]
#     # save and plot the PCs with sorting
#     np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-index_{folder}_{name}.txt',arg_cosine)
#     np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-cosine_{folder}_{name}.txt',cosine_pca[arg_cosine])
#     np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-pca-components_{folder}_{name}.txt',pca_components_sorted)
#     fig_components_sorted = px.imshow(
#         pca_components_sorted,
#         x=columns, y=[f"PC{i}" for i in arg_cosine],
#         color_continuous_scale="RdBu_r", color_continuous_midpoint=0
#     )
#     fig_components_sorted.write_html(f'{Path(f"{out_fig}")}/condition-sorted-pca_components_{folder}_{name}.html')
#     # plot the cosine 
#     plt.figure(figsize=(15,12))
#     plt.barh(np.arange(num_pc),cosine_pca[arg_cosine[::-1]],align='center')
#     plt.yticks(np.arange(num_pc),[f"PC{i}" for i in arg_cosine[::-1]])
#     plt.xlabel(r"$cos\left<condition\ vector,PC\right>$")
#     plt.title(f"{folder}")
#     plt.savefig(f'{Path(f"{out_fig}")}/condition-sorted-cosine_{folder}_{name}.png')
#     plt.close()


#     # Draw projections onto the PCs
#     for i_pc in range(len(pca_components)):
#         base = pca_components[i_pc]
#         df_pca[f"proj{i_pc}"] = df_pca.apply(lambda x:np.dot(base,x.loc[columns]),axis=1)
#     pc2proj = arg_cosine[:3]
#     df_pca_extremes = df_pca.loc[df_pca["condition"].isin(groups)]

#     # 3d projection
#     figproj = plt.figure(figsize=(15,12))
#     ax = figproj.add_subplot(projection="3d")
#     for condi in groups[::-1]:
#         pc_x = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[0]}"],
#         pc_y = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[1]}"],
#         pc_z = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[2]}"],
#         ax.scatter(
#             pc_x, pc_y, pc_z,
#             s=49,alpha=0.2,label=f"{condi}"
#         )
#     ax.set_xlabel(f"proj {pc2proj[0]}")
#     ax.set_ylabel(f"proj {pc2proj[1]}")
#     ax.set_zlabel(f"proj {pc2proj[2]}")
#     ax.xaxis.pane.set_edgecolor('black')
#     ax.yaxis.pane.set_edgecolor('black')
#     ax.zaxis.pane.set_edgecolor('black')
#     ax.xaxis.pane.fill = False
#     ax.yaxis.pane.fill = False
#     ax.zaxis.pane.fill = False
#     ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#     ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#     ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#     ax.legend(loc=(1.04,1.0))
#     figproj.savefig(f'{Path(f"{out_fig}")}/pca_projection3d_{folder}_{name}_pc{"".join([str(p) for p in pc2proj])}.png')
#     plt.close(figproj)

#     # 3d projections, all conditions
#     for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
#         if condi == groups[1]:
#             continue
#         figproj = plt.figure(figsize=(15,12))
#         ax = figproj.add_subplot(projection="3d")
        
#         pc_x = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[0]}"],
#         pc_y = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[1]}"],
#         pc_z = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[2]}"],
#         ax.scatter(
#             pc_x, pc_y, pc_z,
#             edgecolor='white',facecolor=sns.color_palette('tab10')[0],
#             s=49,alpha=0.2,label=f"{groups[-1]}"
#         )
        
#         pc_x = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[0]}"],
#         pc_y = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[1]}"],
#         pc_z = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[2]}"],
#         ax.scatter(
#             pc_x, pc_y, pc_z,
#             edgecolor='white',facecolor=sns.color_palette('tab10')[list_colors[experiment][d]],
#             s=49,alpha=0.2,label=f"{condi}"
#         )
#         ax.set_xlabel(f"proj {pc2proj[0]}")
#         ax.set_ylabel(f"proj {pc2proj[1]}")
#         ax.set_zlabel(f"proj {pc2proj[2]}")
#         ax.xaxis.pane.set_edgecolor('black')
#         ax.yaxis.pane.set_edgecolor('black')
#         ax.zaxis.pane.set_edgecolor('black')
#         ax.xaxis.pane.fill = False
#         ax.yaxis.pane.fill = False
#         ax.zaxis.pane.fill = False
#         ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#         ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#         ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#         ax.legend(loc=(1.04,0.5))
#         figproj.savefig(f'{Path(f"{out_fig}")}/pca_projection3d_{folder}_{name}_condi-{str(condi).replace(".","-")}_pc{"".join([str(p) for p in pc2proj])}.png')
#         plt.close(figproj)
 
#     # 2d projections
#     sns.set_style("whitegrid")
#     for first,second in ((0,1),(0,2),(1,2)):
#         plt.figure(figsize=(15,12))
#         sns_plot = sns.scatterplot(
#             data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[0])],
#             x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
#             color=sns.color_palette("tab10")[1],s=49,alpha=0.5
#         )
#         sns_plot = sns.scatterplot(
#             data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[1])],
#             x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
#             color=sns.color_palette("tab10")[0],s=49,alpha=0.5
#         )
#         sns_plot.figure.savefig(f'{Path(f"{out_fig}")}/pca_projection2d_{folder}_{name}_pc{pc2proj[first]}{pc2proj[second]}.png')
#         plt.close()

#     # 2d projections, all conditions
#     for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
#         if condi == groups[1]:
#             continue
#         for first,second in ((0,1),(0,2),(1,2)):
#             plt.figure(figsize=(15,12))
#             sns_plot = sns.scatterplot(
#                 data=df_pca[df_pca["condition"].eq(groups[1])],
#                 x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
#                 color=sns.color_palette("tab10")[0],s=49,alpha=0.5
#             )
#             sns_plot = sns.scatterplot(
#                 data=df_pca[df_pca["condition"].eq(condi)],
#                 x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
#                 color=sns.color_palette("tab10")[list_colors[experiment][d]],s=49,alpha=0.5
#             )
#             plt.xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[first]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#             plt.ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[second]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
#             sns_plot.figure.savefig(f'{Path(f"{out_fig}")}/pca_projection2d_{folder}_{name}_condi-{str(condi).replace(".","-")}_pc{pc2proj[first]}{pc2proj[second]}.png')
#             plt.close()

#     return None






























#===========================================================



# for org in organelles:
#     print(org)
#     for exp in exp_names:
#         fig, ax = plt.subplots(1, 1,  facecolor='w', edgecolor='k')
#         df_bycell_exp = df_bycell[df_bycell["experiment"] == exp]
#         data_columns = ["mean", "error", "folder", "org", "condition","exp_name"]
#         df = pd.DataFrame(columns=data_columns)
#         for folder in df_bycell_exp["folder"].unique():
#             df_bycell_check=df_bycell_exp[df_bycell_exp["organelle"].eq(org)  & df_bycell_exp["folder"].eq(folder)]
#             for field in df_bycell_check["field"].unique():
#                 df_bycell_check1=df_bycell_check[df_bycell_check["field"].eq(field)]
#                 meanx=df_bycell_check1[f"{var}"].mean()
#                 stdx=df_bycell_check1[f"{var}"].std()
#                 cond2=df_bycell_check1["condition"].values[0]
#                 condition=df_bycell_check1["condition"].values[0]
#                 exp_name1=df_bycell_check1["experiment"].values[0]
#                 print("mean=",meanx,"error",stdx,"folder=",folder,"org=",org,"cond=",df_bycell_check["condition"].values[0])
#                 plt.errorbar(cond2,meanx,fmt="o")
#                 # plt.errorbar(cond2,meanx, yerr=stdx,fmt="o")
#                 plt.text(cond2,meanx,f"{folder}-{field}")
#             ax.set_title(f"{org}-{exp}-{var} ", fontsize=font_size)



for cond in chx_org["condition"].unique():
    print(cond)
    print(chx_org["org"].unique())
    group1=chx_org.loc[chx_org["condition"]==cond,"mean"]
    group2=group_control
    var1=np.var(group1); var2=np.var(group2)
    if var1>var2:
        print(var1/var2)
    else:
        print(var2/var1)
    
    
    tTest=stats.ttest_ind(a=group1, b=group2, equal_var=False)
    print(tTest)



an.get_independent_measurement(organelles[4],exp_names,df_bycell,"mean",with_error_bar=True)

an.get_notch_plot_independent_measure(organelles,exp_names,df_bycell)


       
#########===========================================================###
###########=======================ploting=======================#######
#########===========================================================###
##==========================================##
#####=======================organelle properties vs condition ==========#
###==========================================#

for org in organelles:
    fig, ax = plt.subplots(1, 1)
    data=[]
    df_condition1=[]
    ver_line=-1
    for exp in exp_names:
        ver_line=ver_line+ len(extremes[f"{exp}"])
        print(colored(f"Experiment={exp}:", "red"))
        df_bycell_exp=df_bycell.loc[df_bycell["organelle"].eq(org) & df_bycell["experiment"].eq(exp)]
        for cond in  df_bycell_exp["condition"].unique():
            df_bycell_chx=df_bycell_exp.loc[ df_bycell_exp["condition"].eq(cond)]
            print("org",org,"experiment:",exp,"condition:",cond)
           
            # data[f"{cond}"]=df_bycell_chx["total-fraction"]
            var="mean"
            data.append(df_bycell_chx[f"{var}"].sample(frac=1.0) ) ## choose frac value for randamoly selecting data
            df_condition1.append(cond)     
        palette = sns.color_palette("bright")
        ax= sns.boxplot(data=data,notch=True,palette=palette) #box plot
        # pairs=list(df_condition)
        # annotator=Annotator(ax,pairs,data=data)
        # annotator.apply_and_annotate()
        # ax=sns.stripplot(data = data)
        # ax= sns.violinplot(data=data) #violinplot
        # ax = sns.stripplot(data=data) #violinplot
        ax.set_xticklabels(df_condition1)
        # plt.axvline(x = 0.5, color = 'b', label = 'axvline - full height') #to draw blue line vertically
        plt.axvline(x = ver_line+0.5, color = 'r', label = 'axvline - full height')
        plt.xlabel(f"BFA/CHX concentartion",fontsize=font_size)
        plt.ylabel(f"{org}-{var}",fontsize=font_size)
        plt.title(f"{var}- {org} vs BFA/CHX concentration",fontsize=font_size)
    # plt.savefig(f'{Path(f"{out_fig}")}/Oct20_CHX_concentartion_{org}-{var}.png')




        
    # plt.show()

##==========================================##
#####=======================bootstrap analaysis for confidence level ==========#
###==========================================#

xname=pd.DataFrame()
for org in organelles:
    fig, ax = plt.subplots(1, 1)
    data=[]
    for condition in df_condition:
        # fig, ax = plt.subplots(1, 1)
        df_bycell_chx=df_bycell.loc[df_bycell["organelle"].eq(org) & df_bycell["condition"].eq(condition)]
        var="total-fraction"
        # data.append(df_bycell_chx[f"{var}"].sample(frac=1.0) ) ## choose frac value for randamoly selecting data
        data1=[]     
        data1=df_bycell_chx[f"{var}"].sample(frac=1.0).to_numpy()
        
        bootstrap_statistics = []
        n_bootstrap_samples=1000
        for _ in range(n_bootstrap_samples):
            bootstrap_sample = np.random.choice(data1, size=len(data1), replace=True)
            bootstrap_statistic = np.mean(bootstrap_sample)  # Change this to your statistic of interest
            bootstrap_statistics.append(bootstrap_statistic)

        confidence_interval = np.percentile(bootstrap_statistics, [2.5, 97.5]) #2.5 and 97.5 lower and uper bound of percent
        
        print("95% Confidence Interval for the Mean:", confidence_interval)

        # plt.hist(bootstrap_statistics, bins=30, edgecolor='k', alpha=0.7)
        # plt.xlabel('Bootstrap Statistics')
        # plt.ylabel('Frequency')
        # plt.title(f"Bootstrap Analysis for: Org: {org}-condition:{cond}-para:{var}",fontsize=font_size)
        # # plt.title('Bootstrap Analysis')
        # plt.axvline(np.mean(data1), color='r', linestyle='dashed', linewidth=2, label='Original Mean')
        # plt.axvline(confidence_interval[0], color='g', linestyle='dashed', linewidth=2, label='95% CI Lower Bound')
        # plt.axvline(confidence_interval[1], color='g', linestyle='dashed', linewidth=2, label='95% CI Upper Bound')
        # plt.legend()
        # plt.savefig(f'{Path(f"{out_fig}")}/boot/Bootstrap_plot_CHX_concentartion_{org}-{cond}-{var}.png')
        xname[f'DF_{org}_{cond}'] = pd.DataFrame(bootstrap_statistics)
        data.append(bootstrap_statistics)
    palette = sns.color_palette("bright")
    ax= sns.boxplot(data=data,notch=True,palette=palette) #box plot
    # pairs=list(df_condition)
    # annotator=Annotator(ax,pairs,data=data)
    # annotator.apply_and_annotate()
    # ax=sns.stripplot(data = data)
    # ax= sns.violinplot(data=data) #violinplot
    # ax = sns.stripplot(data=data) #violinplot
    ax.set_xticklabels(df_condition)
    plt.axvline(x = 0.5, color = 'b', label = 'axvline - full height') #to draw blue line vertically
    plt.axvline(x = len(extremes["cycloheximide"])+0.5, color = 'r', label = 'axvline - full height')
    
    plt.xlabel(f"BFA/CHX concentartion",fontsize=font_size)
    plt.ylabel(f"{org}-{var}",fontsize=font_size)
    plt.title(f"{var}- {org} vs BFA/CHX concentration",fontsize=font_size)
        
xname.to_csv(f'{Path(f"{out_fig}")}/{var}_boot.csv')

##==========================================##
#####=======================Corelation-heatmap==========#
###==========================================#


# Correlation coefficient
for exp in extremes:
    
    print(exp)
    xx=df_corrcoef[df_corrcoef["experiment"]==exp]
    print(xx["condition"].unique())
    cond1=str(xx["condition"].iloc[0] )
    np_corrcoef = xx.loc[:,columns].to_numpy()
    corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
  
    fig = px.imshow(
            corrcoef,
            x=columns,y=columns,
            color_continuous_scale = "RdBu_r",range_color=[-1,1],
            # title=str(xx["condition"].iloc[0] )
        )
    
    fig.write_html(f'{Path(f"{out_fig}/")}/1ucorrcoef-_{exp}_all.html')

    for condition in xx["condition"].unique():
        fig1, ax1 = plt.subplots(1, 1)
        np_corrcoef = xx.loc[df_corrcoef['condition']==condition,['effective-length','cell-area','cell-volume',*properties]].to_numpy()
        corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
        ax1 = sns.heatmap(corrcoef,xticklabels=columns, yticklabels=columns, cmap="YlGnBu",
                         
                         annot=True,fmt=".1f",linewidth=.5,vmin=-1,vmax=1,square=True
                         )
        ax1.set_title(f"{exp}_Condition:{condition}")
        fig = px.imshow(
                corrcoef,
                x=columns,y=columns,
                color_continuous_scale="RdBu_r",range_color=[-1,1]
            )
        fig.write_html(f'{Path(f"{out_fig}")}/1ucorrcoef-nocond_{exp}_{str(condition)}_indivual.html')


        # x=f"proj{pc2proj[first]}";y=f"proj{pc2proj[second]}"
        # plt.xlabel(f"{x}")
        # plt.ylabel(f"{y}")
        # plt.scatter( data[f"{x}"],data[f"{y}"],color = 'red')
        # plt.scatter( data1[f"{x}"],data1[f"{y}"],color = 'green')
