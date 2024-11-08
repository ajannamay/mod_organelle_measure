# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:23:49 2023

@author: Anang
"""


from termcolor import colored
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from scipy.stats import sem

# import scipy.stats as stats
# from organelle_measure.data import read_results_pixel
# import plotly.express as px

def plot_back_imgae(data):
    
    organelle=data["organelle"]
    laser="blue"
    drug=data["experiment"]
    condition=data["condition"]
    hour=data["hour"]
    field=data["field"]
    
    
    
    name=f"{organelle}_{laser}_whi5_{drug}_{condition}_time_{hour}_field_{field}.tiff"
    
    
    
    return 



def cal_mean(df_bycell_noise,property1):
    # Initialize an empty list to store the data
    data = []

    # Loop through each unique organelle and time to calculate the average
    for time in df_bycell_noise['hour'].unique():
        row = {'hour': time}
        for org in df_bycell_noise['organelle'].unique():
            org1 = df_bycell_noise[(df_bycell_noise['organelle'] == org) & (df_bycell_noise['hour'] == time)]
            avg = org1[f'{property1}'].mean() if not org1.empty else None
            row[org] = avg
        data.append(row)

    # Create a DataFrame from the collected data
    df_avg = pd.DataFrame(data)
    df_avg=df_avg.drop(['non-organelle'],axis=1)
    df_avg = df_avg.sort_values(by='hour').reset_index(drop=True)
    return df_avg





def group_by_plot(df_stats_org,group_category,path_out,plot_type,catagery):
    organelles = df_stats_org['organelle'].unique()
    time=df_stats_org["hour"].max()

    # Iterate over each organelle and create a plot
    for organelle in organelles:
        fig, ax = plt.subplots(figsize=(10, 6))
        groups = df_stats_org[df_stats_org['organelle'] == organelle].groupby(f"{group_category}")
        
        for folder, group in groups:
            ax.errorbar(
                group["hour"], group["mean_cell_org"], yerr=group["std_cell_org"], 
                marker="o", linestyle="-", label=folder, capsize=5
            )
            
            # ax.errorbar(
            #     group["hour"], group["mean_cell_org"],  
            #     marker="o", linestyle="-", label=folder, capsize=5
            # )
        
        ax.set_title(f'Average Cell Organelle Volume over Time for {organelle}')
        ax.set_xlabel('Hour')
        ax.set_ylabel('Average Cell Organelle Volume')
        ax.legend(title='Folder')
        
        plt.savefig(f"{path_out}/plots/{organelle}_{catagery}_{plot_type}_mean_scatter_{time}.png")
        
    



    
    
font_size=12


out_fig="../../../Analysis/2023_analysis/csv2_figure/jan"



def cell_boundry(img1, mask):    
   # Overlay the boundary on top of the original image
   plt.imshow(img1, cmap='gray')
   plt.contour(mask, levels=[0.5], colors='red', linestyles='dotted')
   plt.axis('off')
   plt.show()

    

def run_pca(data,org_names1):
    
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    import pandas as pd
    
   
    Y1=data.loc[:,["minutes"]].reset_index(drop=True)
    
    Y2=data.loc[:,["experiment"]].reset_index(drop=True)
    
   
    
    pca_data=data.loc[:,org_names1]
    
    labels=pca_data.keys()
    
    NPC=pca_data.shape[1]
    

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pca_data)
    
    # scaled_data=df_normalized
    
    PCA = PCA(n_components=NPC)
    components = PCA.fit_transform(scaled_data )
    coeff=PCA.components_
    
    pc_name=[]
    for i in range(NPC):
        pc_name.append(f"PC{i+1}")
    
    # Plot explained variance ratio
    explained_variance_ratio = PCA.explained_variance_ratio_
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)
    
    cumVar = pd.DataFrame(np.cumsum(PCA.explained_variance_ratio_)*100, 
                          columns=["cumVarPerc"])
    expVar = pd.DataFrame(PCA.explained_variance_ratio_*100, columns=["VarPerc"])
  
    componentsDf = pd.DataFrame(data = components, columns = pc_name).reset_index(drop=True)

    componentsDf=pd.concat([componentsDf,Y1, Y2], axis=1)
    
    coeffDf = pd.DataFrame(data = coeff, columns = labels)
    pcaDf =componentsDf
    
    return coeffDf,pcaDf,NPC,labels
    
    



def sns_plot(data,prop):        
    xdata=data[data["experiment"]=="control"]
    data1=data[data["experiment"]!="control"]
    
    num_rows = 3
    num_cols = 3  # You can adjust this based on your preference

    # Create a figure and a grid of subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 8))

    # Flatten the axs array to iterate over it
    axs = axs.flatten()
    
    # fig.tight_layout(pad=3.0)
    
    
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.4)
    
    ### delete the extra plots
    for j in range(len(data['organelle'].unique()), len(axs)):
        fig.delaxes(axs[j])

    
    for i, org in enumerate(data1['organelle'].unique()):
        ax = axs[i]
        getdata=data1[data1["organelle"]==org]
        datacnt=xdata[xdata["organelle"]==org]
        lines = []
        for exp in getdata["experiment"].unique():
            getdata1=getdata[getdata["experiment"]==exp]
    
            x0=datacnt["growth"].mean()
            y0=datacnt[f"{prop}"].mean()
            x=[]
            y=[]
            
            xx=getdata1["growth"].unique()
            xx.sort()
            for grow in xx:
                getdatagrow=getdata1[getdata1["growth"]==grow]
                x1=getdatagrow["growth"].mean()
                y1_prop=getdatagrow[getdatagrow["organelle"]==org]
                y1=y1_prop[f"{prop}"].mean()
                x.append(x1)
                y.append(y1)
                
            # ax.scatter(x, y)
            x.append(x0)
            y.append(y0)
            yerr=sem(y)
            
            
            
            line=ax.errorbar(x, y, yerr=yerr, fmt='-o')
        
            # ax.gca().invert_yaxis()
            ax.set_xlabel("growth")
            ax.set_ylabel(f"{prop}")
            lines.append(line)
            
            ax.set_title(f"Organelle={org}")
        ax.invert_xaxis()
    fig.legend(lines, data1['experiment'].unique(),bbox_to_anchor=[0.6, 0.2], loc='lower right')
    

x_rapa=[0.5086999325251468,
0.4299795479981529,
0.4207304260796529,
0.2068316673950961,
0.2011211092472894]
y_rapa=[0.04673339369125918,
0.052724647878076764,
0.06625808728214215,
0.11384754291440298,
0.10415039894961424]


# plt.plot(x_rapa,y_rapa)


def reverse_sns_plot(data,prop):        
    xdata=data[data["experiment"]=="control"]
    data1=data[data["experiment"]!="control"]
    
    
    df_chx=pd.DataFrame()
    df_bfa=pd.DataFrame()
    
    num_rows = 3
    num_cols = 3  # You can adjust this based on your preference

    # Create a figure and a grid of subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 8))

    # Flatten the axs array to iterate over it
    axs = axs.flatten()
    
    # fig.tight_layout(pad=3.0)
    
    
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.4)
    
    ### delete the extra plots
    for j in range(len(data['organelle'].unique()), len(axs)):
        fig.delaxes(axs[j])

    
    for i, org in enumerate(data1['organelle'].unique()):
        ax = axs[i]
        getdata=data1[data1["organelle"]==org]
        datacnt=xdata[xdata["organelle"]==org]
        lines = []
        for exp in getdata["experiment"].unique():
            getdata1=getdata[getdata["experiment"]==exp]
    
            x0=datacnt["growth"].mean()
            y0=datacnt[f"{prop}"].mean()
            y01=datacnt[f"{prop}"]
            se=y0/np.sqrt(len(y01))
            x=[x0]
            y=[y0]
            yerr=[se]
            
            xx=getdata1["growth"].unique()
            xx=sorted(xx, reverse=True)
            for grow in xx:
                getdatagrow=getdata1[getdata1["growth"]==grow]
                x1=getdatagrow["growth"].mean()
                y1_prop=getdatagrow[getdatagrow["organelle"]==org]
                y1=y1_prop[f"{prop}"].mean()
                
                y11=y1_prop[f"{prop}"]
                
                x.append(x1)
                y.append(y1)
                se=y1/np.sqrt(len(y11))
                yerr.append(se)
                
            # ax.scatter(x, y)
            # x.append(x0)
            # y.append(y0)
            # yerr=sem(y)
            
            
            
            print(f'={org}_{exp}_{grow}====================')
            print(y)
            
            line=ax.errorbar(x, y, yerr=yerr, fmt='-o')
            if org=='vacuole' and prop=='total-fraction':
                line=ax.errorbar(x_rapa, y_rapa, yerr=0, fmt='-o')
            
            
            
            if exp=='cycloheximide':
                df_chx[f'{exp}_{org}_growth']=x
                df_chx[f'{exp}_{org}_{prop}']=y
                df_chx[f'{exp}_{org}_error']=yerr
            else:
                df_bfa[f'{exp}_{org}_growth']=x
                df_bfa[f'{exp}_{org}_{prop}']=y
                df_bfa[f'{exp}_{org}_error']=yerr
            
          
            
            ax.set_xlabel("growth")
            ax.set_ylabel(f"{prop}")
            lines.append(line)
            
            ax.set_title(f"Organelle={org}")
    
    
    fig.legend(lines, data1['experiment'].unique(),bbox_to_anchor=[0.6, 0.2], loc='lower right')
    return df_bfa,df_chx
   


def plot_cell(data,prop,growth_rate):
    growth=growth_rate["growth_bfa"]

    fig, ax = plt.subplots(1, 1)
    x_values=[]; y_values=[]; error_values=[]
    for condition in sorted(data["condition"].unique()):
        data1=data[data["condition"]==condition]
        xx=data1[f"{prop}"].mean()
        yy=growth[f"{condition}"]
        error=sem(data1[f"{prop}"])
        x_values.append(xx)
        y_values.append(yy)
        error_values.append(error)
    plt.errorbar(y_values,x_values, yerr=error_values,fmt="o")
    plt.plot(y_values,x_values)
    ax.set_title(f"property:{prop}")
    ax.set_xlabel("growth rate")
    ax.set_ylabel(f"{prop}")
    
        
        
        # print("condition:",yy,"mean:",xx,"error:",error)
        # plt.errorbar(yy,xx, yerr=error,fmt="o")
        # plt.plot(yy,xx)
        # ax.set_title(f"property:{prop}")
        # ax.set_xlabel("growth rate")
        # ax.set_ylabel(f"{prop}")



def plot_organelle(data,organelles,prop,growth_rate):
    growth=growth_rate["growth_chx"]
    values=list(data["condition"].unique())
    values.sort()    
    value_to_number = {value: number for number, value in enumerate(values, 1)}
    # np.std(data1[f"{prop}"])/np.sqrt(len(data1[f"{prop}"]))
    for organelle in organelles:
        fig, ax = plt.subplots(1, 1)
        data1=data[data["organelle"]==organelle]
        for condition in sorted(data["condition"].unique()):
            data2=data1[data1["condition"]==condition]
            xx=data2[f"{prop}"].mean()
            yy=growth[f"{condition}"]
            error=sem(data1[f"{prop}"])
            print("condition:",condition,"mean:",xx,"error:",error)
            plt.errorbar(yy,xx, yerr=error,fmt="o")
            ax.set_title(f"organelle:{organelle}---property:{prop}")
            ax.set_xlabel("growth rate")
            ax.set_ylabel(f"{prop}")



              
def custom_function(row):
    extremes = {
        "cycloheximide":                    [ 0.01,0.1,0.5,2],
        "bfa":           [ 25,75, 100.],
        "control":[0.],
        "hu":[4,1.2]
 
    }

    # Example: Square the value from the Reference_Column
    if row['condition'] in extremes["bfa"]:
        x="bfa"
    if row['condition']==0.0:
        x="control"
       
    if row['condition'] in extremes["cycloheximide"]:
        x="cycloheximide"
        
    if row['condition'] in extremes["hu"]:
        x="hu"
    return x

def growth_rate_function(row):
    extremes = {
        "cycloheximide": [ 0.01,0.1,0.5,2],
        "control":[0.],
        "bfa":[25,75,100],
        "hu":[1.2,4]
 
    }

    if row['condition'] in extremes["control"]:
        if row["condition"]==0:
            x=0.5
    if row['condition'] in extremes["cycloheximide"]:
        if row["condition"]==0.01:
            x=0.4
        if row["condition"]==0.1:
            x=0.29        
        if row["condition"]==0.5:
            x=0.16
        if row["condition"]==2:
            x=0.13
    if row['condition'] in extremes["bfa"]:
        if row["condition"]==25:
            x=0.4
        if row["condition"]==75:
            x=0.29
        if row["condition"]==100:
            x=0.13
    if row['condition'] in extremes["hu"]:
        if row["condition"]==1.2:
            x=0.4
        if row["condition"]==4:
            x=0.29
    return x
####===========================

##==========================================##
#####=======================PCA-Analysis==========#
###==========================================#

from pathlib import Path

def make_pca_plots1(df_bycell,organelles,experiment,out_fig,list_colors,property,
                    has_volume=False,is_normalized=True,non_organelle=False):
    # df_bycell=data
    # experiment="cycloheximide"
    # property="total-fraction"
    data1=df_bycell[df_bycell["experiment"]=="control"]
    data2=df_bycell[df_bycell["experiment"]==f"{experiment}"]
    groups=[data1["condition"].max(),data2["condition"].max()]

    folder1 =experiment# experiments[experiment]
    name = f"{'all-conditions' if groups is None else 'extremes'}_\
        {'has-cytoplasm' if non_organelle else 'no-cytoplasm'}_\
            {'cell-volume' if has_volume else 'organelle-only'}_{property}_{'norm-mean-std' if is_normalized else 'raw'}"
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
    
    df_pca.reset_index(inplace=True)

    # Find the the direction of the condition change:
    df_centroid = df_pca.groupby("condition")[list(columns.keys())].mean()
    fitter_centroid = LinearRegression(fit_intercept=False)
    np_centroid = df_centroid.to_numpy()
    fitter_centroid.fit(np_centroid,np.ones(np_centroid.shape[0]))
    
    vec_centroid_start = df_centroid.loc[groups[0],:].to_numpy()
    vec_centroid_start[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_start[:-1]))/fitter_centroid.coef_[-1]

    vec_centroid_ended = df_centroid.loc[groups[-1],:].to_numpy()
    vec_centroid_ended[-1] = (1 - np.dot(fitter_centroid.coef_[:-1],vec_centroid_ended[:-1]))/fitter_centroid.coef_[-1]

    vec_centroid = vec_centroid_ended - vec_centroid_start
    vec_centroid = vec_centroid/np.linalg.norm(vec_centroid)
    np.savetxt(f'{Path(f"{out_fig}")}/condition-vector_{folder1}_{name}.txt',vec_centroid)
    

   # Get Principal Components (PCs)
    np_pca = df_pca[columns].to_numpy()
    pca = PCA(n_components=num_pc)
    pca.fit(np_pca)
    pca_components = pca.components_
    pca_var_ratios = pca.explained_variance_ratio_
    
    # ###Anang pca
    # # Perform PCA
    # X_an=np_pca
    # X_r = pca.fit(X_an).transform(X_an)

    # Percentage of variance explained by each of the selected components
    explained_variance = pca.explained_variance_ratio_


    columns=list(columns.keys())
    import plotly.express as px
    # Calculate cos<condition,PCs>, and realign the PCs
    cosine_pca = np.dot(pca_components,vec_centroid)
    for c in range(len(cosine_pca)):
        if cosine_pca[c] < 0:
            pca_components[c] = -pca_components[c]
    cosine_pca = np.abs(cosine_pca)
    # save and plot PCs without sorting.
    np.savetxt(f'{Path(f"{out_fig}")}/cosine_{folder1}_{name}.txt',cosine_pca)
    np.savetxt(f'{Path(f"{out_fig}")}/pca-components_{folder1}_{name}.txt',pca_components)
    np.savetxt(f'{Path(f"{out_fig}")}/pca-explained-ratio_{folder1}_{name}.txt',pca_var_ratios)
    fig_components = px.imshow(
        pca_components,
        x=columns, y=[f"PC{i}" for i in range(num_pc)],
        color_continuous_scale="RdBu_r", color_continuous_midpoint=0
    )
    fig_components.write_html(f'{Path(f"{out_fig}")}/pca_components_{folder1}_{name}.html')

    # Sort PCs according to the cosine
    arg_cosine = np.argsort(cosine_pca)[::-1]
    pca_components_sorted = pca_components[arg_cosine]
    # save and plot the PCs with sorting
    np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-index_{folder1}_{name}.txt',arg_cosine)
    np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-cosine_{folder1}_{name}.txt',cosine_pca[arg_cosine])
    np.savetxt(f'{Path(f"{out_fig}")}/condition-sorted-pca-components_{folder1}_{name}.txt',pca_components_sorted)
    fig_components_sorted = px.imshow(
        pca_components_sorted,
        x=columns, y=[f"PC{i}" for i in arg_cosine],
        color_continuous_scale="RdBu_r", color_continuous_midpoint=0
    )
    fig_components_sorted.write_html(f'{Path(f"{out_fig}")}/condition-sorted-pca_components_{folder1}_{name}.html')
    # plot the cosine 
    plt.figure(figsize=(15,12))
    plt.barh(np.arange(num_pc),cosine_pca[arg_cosine[::-1]],align='center')
    plt.yticks(np.arange(num_pc),[f"PC{i}" for i in arg_cosine[::-1]])
    plt.xlabel(r"$cos\left<condition\ vector,PC\right>$")
    plt.title(f"{folder1}")
    plt.savefig(f'{Path(f"{out_fig}")}/condition-sorted-cosine_{folder1}_{name}.png')
    plt.close()


    # Draw projections onto the PCs
    for i_pc in range(len(pca_components)):
        base = pca_components[i_pc]
        df_pca[f"proj{i_pc}"] = df_pca.apply(lambda x:np.dot(base,x.loc[columns]),axis=1)
    pc2proj = arg_cosine[:3]
    df_pca_extremes = df_pca.loc[df_pca["condition"].isin(groups)]

    # 3d projection
    figproj = plt.figure(figsize=(15,12))
    ax = figproj.add_subplot(projection="3d")
    for condi in groups[::-1]:
        print(condi)
        pc_x = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[0]}"],
        pc_y = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[1]}"],
        pc_z = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[2]}"],
        
        
        # pc_x = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj{pc2proj[0]}"],
        # pc_y = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj4"],
        # pc_z = df_pca_extremes.loc[df_pca_extremes["condition"].eq(condi),f"proj2"],
        
        
        ax.scatter(
            pc_x, pc_y, pc_z,
            s=5,alpha=0.2,label=f"{condi}"
        )
    ax.set_xlabel(f"proj {pc2proj[0]}")
    ax.set_ylabel(f"proj {pc2proj[1]}")
    ax.set_zlabel(f"proj {pc2proj[2]}")
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
    ax.legend(loc=(1.04,1.0))
    figproj.savefig(f'{Path(f"{out_fig}")}/pca_projection3d_{folder1}_{name}_pc{"".join([str(p) for p in pc2proj])}.png')
    # plt.close(figproj)

    # 3d projections, all conditions
    for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
        print(d,condi)
        if condi == groups[1]:
            continue
        figproj = plt.figure(figsize=(15,12))
        ax = figproj.add_subplot(projection="3d")
        
        pc_x = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[0]}"],
        pc_y = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[1]}"],
        pc_z = df_pca.loc[df_pca["condition"].eq(groups[-1]),f"proj{pc2proj[2]}"],
        ax.scatter(
            pc_x, pc_y, pc_z,
            edgecolor='white',facecolor=sns.color_palette('tab10')[0],
            s=49,alpha=0.2,label=f"{groups[-1]}"
        )
        
        pc_x = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[0]}"],
        pc_y = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[1]}"],
        pc_z = df_pca.loc[df_pca["condition"].eq(condi),f"proj{pc2proj[2]}"],
        ax.scatter(
            pc_x, pc_y, pc_z,
            edgecolor='white',facecolor=sns.color_palette('tab10')[list_colors[experiment][d]],
            s=49,alpha=0.2,label=f"{condi}"
        )
        ax.set_xlabel(f"proj {pc2proj[0]}")
        ax.set_ylabel(f"proj {pc2proj[1]}")
        ax.set_zlabel(f"proj {pc2proj[2]}")
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.zaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.set_xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[0]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.set_ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[1]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.set_zlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[2]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
        ax.legend(loc=(1.04,0.5))
        figproj.savefig(f'{Path(f"{out_fig}")}/pca_projection3d_{folder1}_{name}_condi-{str(condi).replace(".","-")}_pc{"".join([str(p) for p in pc2proj])}.png')
        # plt.close(figproj)
 
    # 2d projections
    sns.set_style("whitegrid")
    for first,second in ((0,1),(0,2),(1,2)):
        plt.figure(figsize=(15,12))
        data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[0])]
        data1=df_pca_extremes[df_pca_extremes["condition"].eq(groups[1])]

        l1=data["condition"].values[0]
        l2=data1["condition"].values[0]
        sns_plot = sns.scatterplot(
            data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[0])],
            x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
            color=sns.color_palette("tab10")[1],s=49,alpha=0.5,
            
        )
        
        sns_plot = sns.scatterplot(
            data=df_pca_extremes[df_pca_extremes["condition"].eq(groups[1])],
            x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
            color=sns.color_palette("tab10")[0],s=49,alpha=0.5
        )
        plt.legend([f"{l1}" , f"{l2}"])
        sns_plot.figure.savefig(f'{Path(f"{out_fig}")}/pca_projection2d_{folder1}_{name}_pc{pc2proj[first]}{pc2proj[second]}.png')
        plt.close()

    # 2d projections, all conditions
    for d,condi in enumerate(np.sort(df_pca["condition"].unique())):
        # print(d,condi)
        if condi == groups[1]:
            print(condi)
            continue
        for first,second in ((0,1),(0,2),(1,2)):

            data=df_pca[df_pca["condition"].eq(groups[1])]
            data1=df_pca[df_pca["condition"].eq(condi)]
            l1=data["condition"].values[0]
            l2=data1["condition"].values[0]
            plt.figure(figsize=(15,12))
            sns_plot = sns.scatterplot(
                data=df_pca[df_pca["condition"].eq(groups[1])],
                x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
                color=sns.color_palette("tab10")[0],s=49,alpha=0.5
            )
            sns_plot = sns.scatterplot(
                data=df_pca[df_pca["condition"].eq(condi)],
                x=f"proj{pc2proj[first]}",y=f"proj{pc2proj[second]}",
                color=sns.color_palette("tab10")[list_colors[experiment][d]],s=49,alpha=0.5
            )
            plt.xlim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[first]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
            plt.ylim(*(np.percentile(df_pca_extremes[f"proj{pc2proj[second]}"].to_numpy(),[1,99])+np.array([-0.1,0.1])))
            plt.legend([f"{l1}" , f"{l2}"])
            sns_plot.figure.savefig(f'{Path(f"{out_fig}")}/pca_projection2d_{folder1}_{name}_condi-{str(condi).replace(".","-")}_pc{pc2proj[first]}{pc2proj[second]}.png')
            plt.close()
            
    return df_pca,explained_variance
   
####get details about each experimets:            
def get_details(exp_names, df_bycell):
    if isinstance(exp_names, list):
        # If exp_names is a list, iterate over each experiment name
        for exp in exp_names:
            process_get_details(exp, df_bycell)
    else:
        # If exp_names is a single value, process it as one experiment
        process_get_details(exp_names, df_bycell)
def process_get_details(exp_name, df_bycell):
    df_bycell_volume = df_bycell[df_bycell["experiment"] == exp_name]
    print(colored(f"Experiment: {exp_name}", "red", attrs=["bold"]))
    for cond in df_bycell_volume["condition"].unique():
        df_bycell_volume1 = df_bycell_volume.loc[df_bycell_volume["condition"].eq(cond)]
        
        xx= df_bycell_volume1["growth"].unique()
        
        unique_folder = df_bycell_volume1["folder"].unique()
        print(
            colored(f"Condition={cond}:", "red"),
            colored(unique_folder, "green"),
            ", Times:",
            colored(len(unique_folder), "red")
        )
        # print(df_bycell_volume1)
        print("Number of cells:", df_bycell_volume1.shape[0] / 7)
        
        print(
            colored(f"Growth={xx}:", "blue"),
           
        )
        
        
        

####cell voume plots            
def cell_volume_notch(extremes,df_bycell):
    if isinstance(extremes, dict):
        fig, ax = plt.subplots(1, 1)
        data=[]
        drug_condition=[]
        # If exp_names is a list, iterate over each experiment name
        for exp in extremes:
            process_cell_volume(exp, df_bycell,ax,data,drug_condition)
    else:
        fig, ax = plt.subplots(1, 1)
        data=[]
        drug_condition=[]
        # If exp_names is a single value, process it as one experiment
        process_cell_volume(extremes, df_bycell,ax,data,drug_condition)
def process_cell_volume(exp, df_bycell,ax,data,drug_condition):
    df_bycell_chx=df_bycell.loc[ df_bycell["experiment"].eq(exp)]
    for condition in df_bycell_chx["condition"].unique():     
        if condition ==0:
            val=np.median(df_bycell_chx["cell-volume"])
            plt.axhline(y = val, color = 'b', label = 'axvline - full height')
            print(f"experiment: {exp}", val)
        data.append(df_bycell_chx["cell-volume"])
        val=np.median(df_bycell_chx["cell-volume"])
        print(f"experiment: {exp}, Condition: {condition}", ",cell volume:", val)
        drug_condition.append(condition)
        ax= sns.boxplot(data=data,notch=True) #box plot  
    ax.set_xticklabels(drug_condition)
    
    plt.xlabel(f"{exp} concentartion",fontsize=font_size)
    plt.ylabel("cell-volume",fontsize=font_size)
    plt.title("cell-volume variation",fontsize=font_size)


####get_independent_measurement          
def get_independent_measurement(organelles,exp_names,df_bycell,var,with_error_bar=False):
    """
    Parameters
    ----------
    organelles : A dictionay of organelles, or as single name.
    exp_names : list of experiments like hu, cycliheximide, bfa.
    df_bycell :data frame of all experiments.
    var: A string the paramenter you want to plot
    Returns
    error bar: Set to True to include error bars (default is False)
    -------
    None.

    """
    if isinstance(organelles, list):
        # If exp_names is a list, iterate over each experiment name
        for org in organelles:
           df= process_get_independent_measurement(org,exp_names,df_bycell,var,with_error_bar)
            
    else:
        
        df=process_get_independent_measurement(organelles,exp_names,df_bycell,var,with_error_bar)
    return df

def process_get_independent_measurement(org,exp_names,df_bycell,var,with_error_bar):
    # org=organelles[3]
    fig, ax = plt.subplots(1, 1,  facecolor='w', edgecolor='k')
    for exp in exp_names:
        data1 = df_bycell[df_bycell["experiment"] == exp]
        data_columns = ["mean", "error", "folder", "org", "condition","exp_name"]
        df = pd.DataFrame(columns=data_columns)
        exp_mean=[]
                
        values=list(df_bycell["condition"].unique())
        values.sort()    
        value_to_number = {value: number for number, value in enumerate(values, 1)}
    
        for folder in data1["folder"].unique(): ##each day experiment
            
            data2=data1[data1["organelle"].eq(org)  & data1["folder"].eq(folder)]
            meanx=data2[f"{var}"].mean()
            stdx=data2[f"{var}"].std()
            condition=data2["condition"].values[0]
            key=value_to_number[condition]
            
            exp_name1=data2["experiment"].values[0]
            

            exp_mean.append(meanx)
            
            # print("measure=",f'{var}',"mean=",meanx,"error",stdx,"folder=",folder,"org=",org,"cond=",df_bycell_check["condition"].values[0])
            
            data = {
        "mean": [meanx],
        "error": [stdx],
        "folder": [folder],
        "org": [org],
        "condition": [condition],
        "exp_name":[exp_name1]
    }
            df = df.append(pd.DataFrame(data), ignore_index=True)
            if with_error_bar:
                # Create a scatter plot with error bars
                plt.errorbar(key,meanx, yerr=stdx,fmt="o")   
            else:
                # Create a regular scatter plot
                plt.scatter(key,meanx)
            plt.text(key+0.05,meanx,f"{condition}",color="red") #{folder}:
 
        ax.set_title(f"{org}-{var} ", fontsize=font_size)
        
        from scipy.stats import sem
        for m_exp in data1["condition"].unique():
            key=value_to_number[m_exp]
            xx=data1[data1["condition"]==m_exp]
            yy=xx[f"{var}"].mean()
            yyer=sem(xx[f"{var}"])
            
            # yyer=xx[f"{var}"].std()
            plt.errorbar(key,yy,yerr=yyer,fmt="*",ecolor="blue")
            plt.scatter(key,yy,marker="*",s=500,alpha=0.5,cmap="blue") #cmap="nipy_spectral
            
            print("condition:",m_exp,"mean:",yy,"st_er:",yyer,"exp:",exp)

            plt.text(key-0.1,yy,f"{m_exp}",color="blue",fontsize=12)
    plt.savefig(f'{Path(f"{out_fig1}")}/{org}-{var}.png')
    
    return df
    
###===================================================================================================================
####get_notch_plot_independent_measure           
def get_notch_plot_independent_measure(organelles,exp_names,df_bycell):
    for org in organelles: ###gpt code
        print("org", org)
        ii = -1
        fig, ax = plt.subplots(1, 4, figsize=(12, 6), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace=0.5, wspace=0.001)
        for exp in exp_names:
            ii = ii + 1
            df_bycell_exp = df_bycell[df_bycell["experiment"] == exp]
            for folder in df_bycell_exp["folder"].unique():
                data = []
                df_condition1 = []
                for cond in df_bycell_exp["condition"].unique():
                    df_bycell_chx = df_bycell_exp.loc[df_bycell_exp["organelle"].eq(org) & df_bycell_exp["condition"].eq(cond) & df_bycell_exp["folder"].eq(folder)]
                    var = "total-fraction"
                    data.append(df_bycell_chx[f"{var}"].sample(frac=1.0))
                    df_condition1.append(cond)
                palette = sns.color_palette("bright")
                ax[ii]=sns.boxplot(data=data, notch=True, palette=palette, ax=ax[ii])
                ax[ii].set_xticklabels(df_condition1)
                ax[ii].set_xlabel(f"{exp}", fontsize=font_size)
                ax[ii].set_title(f" {exp} concentration", fontsize=font_size)
                ax[ii].set_ylabel(f"{org}-{var}", fontsize=font_size)
        plt.tight_layout()
        plt.show()

###=========================================================================================================

####get_each_date_field
            
def get_each_date_field(organelles,exp_names,df_bycell,var):
    for org in organelles:
        data1=df_bycell[df_bycell["organelle"] == org]
        for exp in exp_names:
           data2=data1[data1["experiment"]==exp]
           fig, ax = plt.subplots(1, 1,  facecolor='w', edgecolor='k')
         
           for folder in data2["folder"].unique():
               data3=data2[ data2["folder"].eq(folder)]
               for field in data3["field"].unique():
                   df_bycell_check1=data3[data3["field"].eq(field)]
                   
                   meanx=df_bycell_check1[f"{var}"].mean()
                   stdx=df_bycell_check1[f"{var}"].std()
                   cond2=df_bycell_check1["condition"].values[0]
                   condition=df_bycell_check1["condition"].values[0]
                
                   print("mean=",meanx,"error",stdx,"folder=",folder,"org=",org,"cond=",condition)
                   plt.errorbar(cond2,meanx,fmt="o")
                   # plt.errorbar(cond2,meanx, yerr=stdx,fmt="o")
                   plt.text(cond2,meanx,f"{folder}-{field}")
               ax.set_title(f"{org}-{exp}-{var} ", fontsize=font_size)
               ax.set_ylabel(f"{var}")
               ax.set_xlabel("CHX concentration")
               
######================================================================================

def corelation_plot(extremes,df_corrcoef,columns,properties,cmap=False):
# Correlation coefficient
    for exp in extremes:
        print(exp)
        fig, ax = plt.subplots(1, 1)
        xx=df_corrcoef[df_corrcoef["experiment"]==exp]
        print(xx["condition"].unique())
        cond1=str(xx["condition"].iloc[0] )
        np_corrcoef = xx.loc[:,columns].to_numpy()
        corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
        if cmap:
            cmap=cmap
            
        else:
            cmap="YlGnBu"
        ax = sns.heatmap(corrcoef,xticklabels=columns, yticklabels=columns, cmap=cmap,
                         
                         annot=True,fmt=".1f",linewidth=.5,vmin=-1,vmax=1,square=True
                         )
        ax.set_title(f"{exp}_all")
        for condition in xx["condition"].unique():
            fig1, ax1 = plt.subplots(1, 1)
            np_corrcoef = xx.loc[df_corrcoef['condition']==condition,['effective-length','cell-area','cell-volume',*properties]].to_numpy()
            corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
            ax1 = sns.heatmap(corrcoef,xticklabels=columns, yticklabels=columns, cmap=cmap,
                             
                             annot=True,fmt=".1f",linewidth=.5,vmin=-1,vmax=1,square=True
                             )
            ax1.set_title(f"{exp}_Condition:{condition}")
            
            
            
######################################################



###=============================================================================================            
def get_df_corr_measurement(organelles,exp_names,df_bycell,var,with_error_bar=False):
    """
    Parameters
    ----------
    organelles : A dictionay of organelles, or as single name.
    exp_names : list of experiments like hu, cycliheximide, bfa.
    df_bycell :data frame of all experiments.
    var: A string the paramenter you want to plot
    Returns
    error bar: Set to True to include error bars (default is False)
    -------
    None.

    """
    if isinstance(organelles, list):
        # If exp_names is a list, iterate over each experiment name
        for org in organelles:
           df= process_get_df_corr_measurement(exp_names,df_bycell,var,with_error_bar)
            
    else:
        
        df=process_get_df_corr_measurement(exp_names,df_bycell,var,with_error_bar)
    return df

def process_get_df_corr_measurement(exp_names,df_bycell,var,with_error_bar):
    # org=organelles[3]
    fig, ax = plt.subplots(1, 1,  facecolor='w', edgecolor='k')
    for exp in exp_names:
        df_bycell_exp = df_bycell[df_bycell["experiment"] == exp]
        data_columns = ["mean", "error", "folder", "org", "condition","exp_name"]
        df = pd.DataFrame(columns=data_columns)
        exp_mean=[]
                
        values=list(df_bycell["condition"].unique())
        values.sort()    
        value_to_number = {value: number for number, value in enumerate(values, 1)}
    
        for folder in df_bycell_exp["folder"].unique(): ##each day experiment
            
            df_bycell_check=df_bycell_exp[ df_bycell_exp["folder"].eq(folder)]
            meanx=df_bycell_check[f"{var}"].mean()
            stdx=df_bycell_check[f"{var}"].std()
            condition=df_bycell_check["condition"].values[0]
            key=value_to_number[condition]
            
            exp_name1=df_bycell_check["experiment"].values[0]
            

            exp_mean.append(meanx)
            
            print("measure=",f'{var}',"mean=",meanx,"error",stdx,"folder=",folder,"cond=",df_bycell_check["condition"].values[0])
            
            data = {
        "mean": [meanx],
        "error": [stdx],
        "folder": [folder],
        "condition": [condition],
        "exp_name":[exp_name1]
    }
            df = df.append(pd.DataFrame(data), ignore_index=True)
            if with_error_bar:
                # Create a scatter plot with error bars
                plt.errorbar(key,meanx, yerr=stdx,fmt="o")   
            else:
                # Create a regular scatter plot
                plt.scatter(key,meanx)
            plt.text(key+0.05,meanx,f"{folder}:{condition}",color="red")
 
        ax.set_title(f"{exp}-{var} ", fontsize=font_size)
        for m_exp in df["condition"].unique():
            key=value_to_number[m_exp]
            xx=df[df["condition"]==m_exp]
            yy=xx["mean"].mean()
            yyer=xx["mean"].std()
            plt.errorbar(key,yy,yerr=yyer,fmt="*",ecolor="blue")
            plt.scatter(key,yy,marker="*",s=500,alpha=0.5,cmap="blue") #cmap="nipy_spectral
            

            plt.text(key-0.1,yy,f"{m_exp}",color="blue",fontsize=12)
    
    return df



####====count-vs size or total-fraction==========================================================




import pandas as pd

def divide_numeric_columns(df, divisor_column):
    # Select only numeric columns
    numeric_columns = df.select_dtypes(include=['number']).columns
    
    # Create a copy of the original DataFrame
    result_df = df.copy()

    # Iterate over numeric columns and perform division
    for col in numeric_columns:
        try:
            # Perform division
            result_df[col] = df[col] / df[divisor_column]
        except ZeroDivisionError:
            # Handle division by zero exception
            print(f"Warning: Division by zero in column '{col}'")
            result_df[col] = float('nan')  # Set result to NaN in case of division by zero
        except Exception as e:
            # Handle other exceptions
            print(f"Error: {e} in column '{col}'")

    return result_df




###=============================================================================================            
def get_count_measurement(organelles,exp_names,df_bycell,var,with_error_bar=False):
    """
    Parameters
    ----------
    organelles : A dictionay of organelles, or as single name.
    exp_names : list of experiments like hu, cycliheximide, bfa.
    df_bycell :data frame of all experiments.
    var: A string the paramenter you want to plot
    Returns
    error bar: Set to True to include error bars (default is False)
    -------
    None.

    """

    for org in organelles:
        fig, ax = plt.subplots(1, 1,  facecolor='w', edgecolor='k')
        for cond in df_bycell["condition"].unique():
            zz=df_bycell[df_bycell["condition"]==cond]
            xx=zz["count"]
            yy=zz[f"{var}"]
            ax.scatter(xx, yy,label=cond)
        
            # plt.scatter(xx,yy) #cmap="nipy_spectral
            ax.set_ylabel(f"{org}-{var}", fontsize=font_size)
            ax.set_xlabel(f"{org}-count", fontsize=font_size)
    ax.set_title(f"{org}-{var} ", fontsize=font_size)
    plt.show()

            






















         

# PLOTS from Shimon
from plotly.subplots import make_subplots
import plotly.graph_objects as go
def plot_histo_violin(df,prop_y,prop_x="cell_area",savepath="histo-violin.html"):
    df["plotly"] = df[prop_y].rank(pct=True) # `pct` works but fuck why
    fig = make_subplots(rows=1,cols=1)
    fig.append_trace(
        go.Violin(
            y=df["plotly"],
            # y=np.arange(len(df[prop_y].unique())),
            x=df[prop_x]
        ),
        row=1, col=1
    )
        
    wid = 2./len(df["plotly"].unique())
    
    fig.update_traces(orientation='h', side='positive', width=wid, points="all",pointpos=0.5,jitter=0.4,line_color="deepskyblue",marker_color="tomato",marker_opacity=0.1,row=1, col=1)
    # print(df[prop_y].unique())
    # print(df["plotly"].unique())
    fig.update_layout(
                    title=Path(savepath).stem.replace("_"," "),
                    yaxis_title=prop_y,
                    yaxis_ticktext=df[prop_y].unique(),
                    yaxis_tickvals=df["plotly"].unique(),
                    xaxis_title=f'{prop_x}'
    )
    fig.show()
    fig.write_html(str(savepath))
    return None

def plot_histo_box(df,prop_y,prop_x="cell_area",savepath=""):
    fig = make_subplots(rows=1,cols=1)
    for prop in df[prop_y].unique():
        fig.append_trace(
            go.Box(x=df.loc[df[prop_y]==prop,prop_x],
                   boxmean='sd'),
            row=1, col=1
        )
    
    fig.update_layout(
                    title=Path(savepath).stem.replace("_"," "),
                    yaxis_title=prop_y,
                    yaxis_ticktext=df[prop_y].unique(),
                    yaxis_tickvals=df[prop_y].unique(),
                    xaxis_title=f'{prop_x}'
    )
    fig.show()
    fig.write_html(str(savepath))
    return None

def plot_group_hist(df,prop_y,prop_x="cell_area"):
    fig = make_subplots(rows=1,cols=1)
   
    for i,prop in enumerate(df[prop_y].unique()):
        # print(prop)
        print(i)
        fig.append_trace(
            go.Histogram(
                x=df.loc[df[prop_y].eq(prop),prop_x],
                name=df[prop_y].unique()[i]),
            row=1, col=1
        )
    
    fig.update_layout(
                    barmode='overlay',
                    
                    yaxis_title=prop_y,
                    xaxis_title=f'{prop_x}',
                    bargap=0.001,bargroupgap=0.05
    )
    fig.update_traces(opacity=0.75)
    # fig.show()
 
    return None


def an_hist_plot(organelles,exp_names,df,var):
    for org in organelles:
        fig, ax = plt.subplots(1, 1)
        for exp in exp_names:
            x=df[f"{var}"]
            for condition in df["condition"].unique():
                
                n_bins=10
                ax.hist(x, bins = n_bins)
                
            
   
import pandas as pd

# Create a DataFrame
df = pd.DataFrame({'Category': ['A', 'B', 'A', 'B'], 'Value': [1, 2, 3, 4]})

# Use groupby and sum
result = df.groupby('Category')['Value'].sum().reset_index()

