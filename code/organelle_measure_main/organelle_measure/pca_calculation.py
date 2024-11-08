# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 10:13:58 2023

@author: Anang
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import matplotlib.colors as colors
###############==========================



        
def pca_viz(pca, X, X_features, sort=False, signed=True):

    pca.fit(X)
    n_rows, n_cols = len(X_features), len(pca.components_)

    inv = pca.inverse_transform(np.eye(n_cols))
    data = pd.DataFrame(inv, columns=X_features).transpose()

    # make zero-values just very, very small
    epsilon = 1e-16
    num_zeros = (data == 0).sum().sum()
    if num_zeros:
        data.loc[:, :] = np.where(data == 0, 1e-16, data)

    # sort by principal component contribution
    if sort:
        data['sum'] = data.abs().sum(axis=1)
        data.sort_values(by='sum', ascending=False, inplace=True)
        data.drop('sum', axis=1, inplace=True)

    # configure heatmap and log-scale colorbar
    scale_min = max(data.abs().min().min(), 1e-2)
    scale_max = data.abs().max().max()
    if signed:
        cmap = 'RdBu_r'
        scaler = colors.SymLogNorm(vmin=-scale_max, vmax=scale_max,
                                   linscale=1.0, linthresh=scale_min)
    else:
        data = data.abs()
        cmap='Reds'
        scaler = colors.LogNorm(vmin=scale_min, vmax=scale_max)

    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(data, cmap=cmap)
    im.set_norm(scaler)
    fig.colorbar(im)

    # configure ticks and labels
    ylabels = data.index
    ax.set_yticklabels(ylabels, fontsize=6, rotation=0, va='center')
    ax.set_yticks(np.arange(n_rows))
    ax.set_ylabel('Original Features')
    xlabels = np.arange(n_cols)
    ax.set_xticklabels(xlabels, fontsize=6, rotation=90, ha='center')
    ax.set_xticks(np.arange(n_cols))
    ax.set_xlabel('Principal Components')

    return fig, ax    

    



def pca_plot(data,pc1,pc2,pc3,Npc,org_names,property1="total-fraction"):
        
    
    # property1="total-fraction"
    
    df_pca = data.pivot_table(index=["experiment", "condition",
                                     "field", "idx-cell", "folder","cell-volume","growth"],
                              columns="organelle", values=f"{property1}", aggfunc="first").reset_index()
    
    # df_pca=df_pca.dropna()
    # exp=df_pca["experiment"].unique()
    # gr=df_pca["growth"].unique()
    
    X=df_pca.loc[:,[*org_names]]
    Y=df_pca.loc[:,["experiment"]]
    

    scaler = StandardScaler()
    
    X = scaler.fit_transform(X)
    
    from sklearn.decomposition import PCA
    PCA = PCA(n_components=Npc)
    components = PCA.fit_transform(X)
    PCA.components_
    
    columns=[]
    for i in range(Npc):
        columns.append(f"PC{i+1}")
    
    # Plot explained variance ratio
    explained_variance_ratio = PCA.explained_variance_ratio_
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)
    
    cumVar = pd.DataFrame(np.cumsum(PCA.explained_variance_ratio_)*100, 
                          columns=["cumVarPerc"])
    expVar = pd.DataFrame(PCA.explained_variance_ratio_*100, columns=["VarPerc"])
    
        
    componentsDf = pd.DataFrame(data = components, columns = columns)
    
    Y=df_pca.loc[:,["growth"]]
    Y1=df_pca.loc[:,["experiment"]]
    pcaDf = pd.concat([componentsDf, Y,Y1], axis=1)
    
    
    def biplot(score,coeff,labels=None):
        xs = score[:,0]
        ys = score[:,1]
        n = coeff.shape[0]
        scalex = 1.0/(xs.max() - xs.min())
        scaley = 1.0/(ys.max() - ys.min())
        plt.scatter(xs * scalex,ys * scaley,s=5)
               
        
        for i in range(n):
            plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
            if labels is None:
                plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'green', ha = 'center', va = 'center')
            else:
                plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')
     
        plt.xlabel("PC{}".format(1))
        plt.ylabel("PC{}".format(2))
        plt.grid()
        
    plt.figure(figsize=(12, 6))
    biplot(components, np.transpose(PCA.components_), list(org_names))
    plt.title("Biplot")
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, marker='o')
    plt.xlabel('Number of Principal Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.title('Cumulative Explained Variance vs. Number of Principal Components')
    plt.grid(True)
    plt.show()
    
    ########################2D PC plot
    
    plt.figure(figsize=(12,10))
    with sns.plotting_context("talk",font_scale=1.25):
        sns.scatterplot(x=f"PC{pc1}", y=f"PC{pc2}",
                        data=pcaDf , 
                        hue="experiment",
                        style="growth",
                        s=10)
        plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
        plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")
        plt.title("PCA with Scaled Penguins Data")
     
        
    plt.figure() 
    pcamodel=PCA    
    ax = sns.heatmap(pcamodel.components_,
                 cmap='YlGnBu',
                 yticklabels=[ "PCA"+str(x) for x in range(1,pcamodel.n_components_+1)],
                 xticklabels=list(org_names),
                 cbar_kws={"orientation": "vertical"})
    ax.set_aspect("equal")
    ax.set_title("Effect of variables on each components")
    
    import plotly.io as pio
    import plotly.express as px
    pio.renderers.default='browser'
    

    
    # fig = px.scatter_3d(pcaDf, x=f"PC{pc1}", y=f"PC{pc2}", z=f"PC{pc3}",
    #           color="growth",symbol='experiment')
    
    fig = px.scatter_3d(pcaDf, x=f"PC{pc1}", y=f"PC{pc2}", z=f"PC{pc3}",
              color='experiment')
    fig.update_traces(marker_size = 3)
       
    fig.show()
    
    
    from sklearn.decomposition import PCA
    
    ##{https://www.brentaustgen.com/blogs/pca-visualization/}
    
    
    X=df_pca.loc[:,[*org_names]]
    X_features=org_names
    
    
    pca = PCA()
    fig, ax = pca_viz(pca, X, X_features, signed=True)
    plt.subplots_adjust(left=0.25, right=0.90, bottom=0.20, top=0.80)
    ax.set_title(" original data")
    plt.show()
    
    
    X_ss = StandardScaler().fit_transform(X)
    fig, ax = pca_viz(pca, X_ss, X_features, signed=True)
    ax.set_title(" scaled data")
    plt.subplots_adjust(left=0.25, right=0.90, bottom=0.20, top=0.80)
    plt.show()
    


        
    ##########=================3D plot
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # pcaDf['experiment']=pd.Categorical(pcaDf['experiment'])
    # my_color=pcaDf['experiment'].cat.codes
    
    # xs = pcaDf[f"PC{pc1}"]
    # ys = pcaDf[f"PC{pc2}"]
    # zs = pcaDf[f"PC{pc3}"]
    # ax.scatter(xs, ys, zs,c=my_color, cmap="Set2_r")
    
    
# from sklearn.decomposition import PCA
# import pandas as pd
# import numpy as np 


# df = df_pca

# array = df.values


# pca_data=df_pca.loc[:,org_names]

# labels=pca_data.keys()

# NPC=pca_data.shape[1]



# Y=df_pca["experiment"]



# #Normalize Data
# def normalize(df):
#     result = df.copy()
#     for feature_name in df.columns:
#         max_value = df[feature_name].max()
#         min_value = df[feature_name].min()
#         result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
#     return result



# df_normalized = normalize(pca_data)

# pca = PCA(n_components = 16)
# pca.fit_transform(df_normalized)
# PCAdf = pd.DataFrame(pca.components_, columns = df_normalized.columns, index = ['PC-1','PC-2','PC-3','PC-4','PC-5','PC-6','PC-7','PC-8','PC-9','PC-10','PC-11','PC-12','PC-13','PC-14','PC-15','PC-16'])
# PCAarray = PCAdf.values

# from sklearn.preprocessing import LabelEncoder
# #Convert all of the "M" class labels as 1, and "B" Labels as 0
# encoder = LabelEncoder()
# encoder.fit(Y)
# encoded_Y = encoder.transform(Y)
# df_v_y_encoded = encoder.transform(df_v_y)

# from sklearn import svm
# #Train again, this time using features from principal component analysis. 
# classifierPCAfeatures = svm.SVC(gamma = "auto", C = 1, kernel = "rbf", decision_function_shape='ovo')
# classifierPCAfeatures = pca.fit(PCAdf, encoded_Y)
# print(classifierPCAfeatures.score(df_v_x, df_v_y_encoded))














    
    
  


'''

    
    pca = PCA().fit(X)
data = pca.inverse_transform(np.eye(X.shape[-1]))
plt.imshow(np.log(data), cmap='hot')
plt.ylabel('Principal Components')
plt.xlabel('Original Features')
plt.colorbar()
plt.show()
    




.

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
            





'''

