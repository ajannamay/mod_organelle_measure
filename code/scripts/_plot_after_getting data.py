# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 16:00:32 2024

@author: Anang
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from organelle_measure.data import read_results
import organelle_measure.vars_allround1data as var
from termcolor import colored
import organelle_measure.an_functions as an
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler #(data-data.mean)/data.std()
import organelle_measure.pca_calculation as pca_plot
import plotly.io as pio
import plotly.express as px
from organelle_measure.tools import load_nd2_plane
from organelle_measure.yeaz import yeaz_preprocesses
pio.renderers.default='browser'

#=================PCA==========================================================




import matplotlib.pyplot as plt
import seaborn as sns



folder=df_org['folder'].unique()

data1=df_org.dropna()

time=df_org['hour'].max()



data1=data[data['folder']=='July01']


property1="total-fraction"
df_pca = data1.pivot_table(index=[ 
                                   "hour","cell-volume",'folder'],
                           columns="organelle", values=f"{property1}", aggfunc="first").reset_index()


org_names1=['peroxisome','ER', 'vacuole', 'golgi', 'mitochondria', 'LD']
pca_data=df_pca.loc[:,org_names1]

pca_data['min']=df_pca['hour']



coeffDf,pcaDf=run_pca(pca_data)

###
pc1=1
pc2=3



plt.figure(figsize=(10, 6))
sns.scatterplot(x=f"PC{pc1}", y=f"PC{pc2}",
                data=pcaDf , 
                # col="experiment",
                style="min",
                hue="min",
                # kind="scatter",
                palette="Set1",
                
                s=30)
plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")

sns.scatterplot(x=f"PC{pc1}", y=f"PC{pc2}",
                data=pcaDf , 
                # col="experiment",
                style="min",
                hue="min",
                # kind="scatter",
                palette="Set1",
                
                s=30)
plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")





# Example DataFrames
data1 = {'x': [1, 2, 3, 4, 5], 'y': [2, 3, 5, 7, 11]}
data2 = {'x': [1, 2, 3, 4, 5], 'y': [1, 4, 9, 16, 25]}

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)



plt.figure(figsize=(10, 6))

# Plot first dataframe
sns.scatterplot(data=df1, x='x', y='y', s=100, label='Dataset 1')

# Plot second dataframe
sns.scatterplot(data=df2, x='x', y='y', s=100, label='Dataset 2')

# Customizing the plot
plt.title('Overlapping Scatter Plot of Two DataFrames')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend(title='Dataset')
plt.show()










   
with sns.plotting_context("paper",font_scale=1.25):
    ax=sns.relplot(x=f"PC{pc1}", y=f"PC{pc2}",
                    data=pcaDf , 
                    # col="experiment",
                    style="min",
                    hue="min",
                    kind="scatter",
                    palette="Set1",
                    
                    s=30)
    plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
    plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")
ax.fig.suptitle(f' Plot for {property1}_{folder}__{time}')
plt.savefig(f"{var.path_out}/plots/pca_{property1}_{folder}_{time}.png")





pallete1=sns.color_palette("bright", 14)




 











 
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(cumulative_explained_variance) + 1), cumulative_explained_variance, marker='o')
# plt.xlabel('Number of Principal Components')
# plt.ylabel('Cumulative Explained Variance')
# plt.title(f'Cumulative Explained Variance vs. Number of Principal Components {property1}')
# plt.grid(True)
# plt.show()

# ###
# pc1=1
# pc2=2
   
# with sns.plotting_context("talk",font_scale=1.25):
#     ax=sns.relplot(x=f"PC{pc1}", y=f"PC{pc2}",
#                     data=pcaDf , 
#                     # col="experiment",
#                     style="hour",
#                     hue="hour",
#                     kind="scatter",
#                     palette="Set1",
                    
#                     s=30)
#     plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
#     plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")
# ax.fig.suptitle(f' Plot for {property1}')
   
# #########################################################################
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler


# property1="total-fraction"
# df_pca = data1.pivot_table(index=[ 
#                                    "hour","cell-volume"],
#                            columns="organelle", values=f"{property1}", aggfunc="first").reset_index()



# df_pca1=df_pca.drop('cell-volume',axis=1)


# data =df_pca1 
# # Assuming the last column is the target
# X = data.drop('hour',axis=1).values
# y = data["hour"].values

# target=y

# # Standardize the data
# scaler = StandardScaler()
# data_scaled = scaler.fit_transform(X)

# # Perform PCA
# pca = PCA(n_components=2)
# principal_components = pca.fit_transform(data_scaled)

# # Create a DataFrame with the principal components
# pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
# pca_df['target'] = y


# centroids = pca_df.groupby('target').mean().reset_index()


# plt.scatter(centroids['PC1'],centroids['PC2'] )
# plt.show()


# centroids.to_csv('centroids.csv')


# with sns.plotting_context("talk", font_scale=1.25):
#     ax = sns.relplot(x=f"PC{pc1}", y=f"PC{pc2}",
#                      data=centroids, 
#                      hue="target",
#                      kind="scatter",
#                       palette="Paired",
#                      s=100)











# target_name=pca_df['target'].unique()

# palette = sns.color_palette("Set1", len(target_name))


# # Plot with Seaborn and Add Centroids
# pc1 = 1
# pc2 = 2

# with sns.plotting_context("talk", font_scale=1.25):
#     ax = sns.relplot(x=f"PC{pc1}", y=f"PC{pc2}",
#                      data=pca_df, 
#                      hue="target",
#                      kind="scatter",
#                      palette="Set1",
#                      s=30)
    
#     for i, target_name in enumerate(target_name):
#         centroid_color = palette[i]
#         plt.scatter(centroids.loc[i, 'PC1'], centroids.loc[i, 'PC2'], 
#                     marker='X', s=100, color=centroid_color, label=f'Centroid {target_name}')
    
#     # for i in range(len(centroids)):
#     #     plt.scatter(centroids.iloc[i, 1], centroids.iloc[i, 2], 
#     #                 marker='X', s=100, color='black')
    
#     plt.xlabel(f"PC{pc1}: " )
#     plt.ylabel(f"PC{pc2}: " )
#     plt.title(f'PCA Plot for Total-fraction')

# plt.show()

# pc1=1
# pc2=2
   
# with sns.plotting_context("talk",font_scale=1.25):
#     ax=sns.relplot(x=f"PC{pc1}", y=f"PC{pc2}",
#                     data=pcaDf , 
#                     # col="experiment",
#                     style="hour",
#                     hue="hour",
#                     kind="scatter",
#                     palette="Set1",
                    
#                     s=30)
#     plt.xlabel(f"PC{pc1}: "+f'{cumulative_explained_variance[int(pc1)]*100:.2f}'+"%")
#     plt.ylabel(f"PC{pc2}: "+f'{cumulative_explained_variance[int(pc2)]*100:.2f}'+"%")
# ax.fig.suptitle(f' Plot for {property1}')


# ax = sns.scatterplot(means[:, 0], means[:, 1], hue=range(4), palette=colors, s=20, ec='black', legend=False, ax=ax)
# plt.show()









# # Plot the results
# plt.figure(figsize=(10, 7))
# colors = ['r', 'g', 'b']
# for target, color in zip(np.unique(target), colors):
#     indices_to_keep = pca_df['target'] == target
#     plt.scatter(pca_df.loc[indices_to_keep, 'PC1'],
#                 pca_df.loc[indices_to_keep, 'PC2'],
#                 c=color,
#                 s=50)
# plt.xlabel('Principal Component 1')
# plt.ylabel('Principal Component 2')
# plt.title('PCA of Iris Dataset')
# # plt.legend(target_names)
# plt.grid()

# plt.show()


# plt.scatter(pca_df[ 'PC1'],
#             pca_df['PC2'],)




#==============================================================================



data=df_bycell_concated



exp_name=["bfa","control","cycloheximide"]
an.get_details(exp_name,data)

data=data[data["cell-volume"]<700]

data["experiment"].unique()




var='total-fraction'

##########histogram plot




for gr in data['experiment'].unique():
  
    data1=data[data['experiment']==gr] 
    
    
    sns.distplot(data1[f'{var}'], hist = True, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 3}, 
                  label = gr)
  
   





for gr in data['experiment'].unique():
  
    data1=data[data['experiment']==gr]  
  
    data1[f'{var}'].hist(bins=50, by=data1['growth'], 
                     alpha=0.8, figsize=(8,8),density=True)



####### sactter plot

for gr in data['organelle'].unique():
  
    data1=data[data['organelle']==gr]  
  
    data1.plot.scatter(x='cell-volume',
                          y='total-fraction',
                          c='condition',
                          colormap='viridis')



fig = px.strip(data, x='experiment', y='total-fraction', color="growth")

# Customize traces
fig.update_traces(jitter=1, opacity=0.8, marker=dict(size=10, line=dict(width=2, color='black')))

fig.show()


###############





import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
                  

fig = px.strip(data, x='experiment', y='total-fraction', color="growth").update_traces(jitter = 1,
                                                                       opacity=0.8,
                                                                       marker_size=10,
                                                                       marker_line_width=1)



# Group and calculate the mean and sem for distance moved
mean = data.groupby('experiment').mean()
sem = data.groupby('experiment').sem()

# Group and calculate the mean and sem for distance moved during dark and light period
mean_period=data.groupby(['experiment','growth']).mean()
sem_period=data.groupby(['experiment','growth']).sem()


# Add traces for mean and sem
fig.add_trace(
    go.Scatter(
        mode='markers',
        x=mean.index, y=data['total-fraction'],
        error_y_array=sem_period['Distance moved'],
        marker=dict(symbol='141', color='rgba(0,0,0,0.6)', size=30,
        line=dict(width=2)
        ),
        showlegend=False
    )
)


#  Customization of y-axis
#fig.update_yaxes(range=[0, 10])

# Figure layout
fig.update_layout(template='simple_white',  width=1000, height=500, title='', yaxis_title='Distance moved [mm]',
                  legend=dict(title='', itemclick='toggle', itemsizing='constant', traceorder='normal',
                  bgcolor='rgba(0,0,0,0)', x=1),
                  #margin=dict(color="black",width=3),
                  xaxis=dict(title='', showticklabels=True, ticks='outside', type='category')
                 )

# Make figure zoomable
config = dict({'scrollZoom':True})

print(df.groupby(['genotype', 'period']).mean())

fig.show(config=config)
























for g in pd.unique(data['experiment']):
    data.loc[data['experiment']==g,'total-fraction'].hist(bins=100,alpha=0.65,by=data['growth'],
                                             label=g,density=True)

    
    # plt.legend(loc='upper left')





for g in pd.unique(data['experiment']):
    data.loc[data['experiment']==g,'total-fraction'].hist(bins=100,alpha=0.65,
                                             label=g,density=True)
plt.legend(loc='upper left')










var='total-fraction'

for org in org_names:
    data1=data[data['organelle']==org]
    fig, ax = plt.subplots(1, 1)
    for exp in exp_names:
        x=data1[data1["experiment"]==exp]
        for condition in x["growth"].unique():
            
            y=x[x['growth']==condition]
            
            z=y[f'{var}']
            
            n_bins=5
            ax.hist(z, bins = n_bins)
            
            
            
            









props=["mean","count","total","cell-volume","total-fraction"]
for prop in props:
   df_bfa,df_chx= an.reverse_sns_plot(data,prop)
   
   
   
########################

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

df = df_bycell_concated


# Data Extraction
# Group and calculate the mean and sem for distance moved during dark and light period
mean_period=df.groupby(['experiment','organelle']).mean()
sem_period=df.groupby(['experiment','organelle']).sem()

# Extract mean from the distance moved during the dark period
dark_mean_A=mean_period['Distance moved'].A['Dark']
dark_mean_B=mean_period['Distance moved'].B['Dark']
dark_mean_C=mean_period['Distance moved'].C['Dark']

# Extract mean from the distance moved during the light period
light_mean_A=mean_period['Distance moved'].A['Light']
light_mean_B=mean_period['Distance moved'].B['Light']
light_mean_C=mean_period['Distance moved'].C['Light']


# Extract sem from the distance moved during the dark period
dark_sem_A=sem_period['Distance moved'].A['Dark']
dark_sem_B=sem_period['Distance moved'].B['Dark']
dark_sem_C=sem_period['Distance moved'].C['Dark']

# Extract sem from the distance moved during the light period
light_sem_A=sem_period['Distance moved'].A['Light']
light_sem_B=sem_period['Distance moved'].B['Light']
light_sem_C=sem_period['Distance moved'].C['Light']


# Scatter plot
fig = px.strip(df, x='genotype', y='Distance moved', color="period",
               color_discrete_sequence=['rgba(0,0,0,0.4)', 'rgba(255,255,0,0.4)']).update_traces(
                                                                       jitter = 1,
                                                                       opacity=0.8,
                                                                       marker_size=10,
                                                                       marker_line_width=1,
                                                                       marker_line_color='rgba(0,0,0,0.8)',
                                                                       #marker_color='rgba(0,0,0,0.8)',
                                                                       showlegend=False)

# Bar graphs with error bars
fig.add_bar(
    name='Dark',
    marker_color='rgba(0,0,0,0.5)', marker_line_color='rgba(0,0,0,0.8)', marker_line_width=1, opacity=0.8,
    x=['A', 'B', 'C'],
    y=[dark_mean_A, dark_mean_B, dark_mean_C],
    error_y=dict(type='data', array=[dark_sem_A, dark_sem_B, dark_sem_C],
                color='rgba(0,0,0,1)', thickness=1.5, width=10)
)

fig.add_bar(
    name='Light',
    marker_color='rgba(255,255,0,0.5)', marker_line_color='rgba(0,0,0,0.8)', marker_line_width=1, opacity=0.8,
    x=['A', 'B', 'C'],
    y=[light_mean_A, light_mean_B, light_mean_C],    
    error_y=dict(type='data', array=[light_sem_A, light_sem_B, light_sem_C],
                color='rgba(0,0,0,1)', thickness=1.5, width=10)
)

# Sample numbers
# Add n numbers
fig.add_trace(go.Scatter(
    x=['B', 'A', 'C'],
    y=[100, 100, 100],
    mode="text",
    name="n numbers",
    text=['n = 30', 'n = 32', 'n = 32'],
    textposition="top center",
    textfont=dict(color='rgba(0,0,0,1)', size=13),
    hoverlabel=dict(bgcolor='white'),
    showlegend=False
))

# Brackets for p-values
# Dark bar 1 to dark bar 2 p-value bracket
x_coords = [0.10, 0.10, 0.428, 0.428]
y_coords = [(dark_mean_A+dark_sem_A)+100, (dark_mean_B+dark_sem_B)+300,
            (dark_mean_B+dark_sem_B)+300, (dark_mean_B+dark_sem_B)+100]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1
    )
    
# Dark bar 1 to dark bar 3 p-value bracket
x_coords = [0.10, 0.10, 0.7674, 0.7674]
y_coords = [(dark_mean_A+dark_sem_A)+700, (dark_mean_A+dark_sem_A)+2000,
            (dark_mean_A+dark_sem_A)+2000, (dark_mean_A+dark_sem_A)+1800]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1
    )
    
# Dark bar 2 to bar 3 p-value bracket
x_coords = [0.437, 0.437, 0.7674, 0.7674]
y_coords = [(dark_mean_B+dark_sem_B)+100, (dark_mean_B+dark_sem_B)+300,
            (dark_mean_B+dark_sem_B)+300, (dark_mean_B+dark_sem_C)+100]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1
    )
    
    
# Light bar 1 to light bar 2 p-value bracket
x_coords = [0.233, 0.233, 0.560, 0.560]
y_coords = [(light_mean_B+light_sem_B)+2500, (light_mean_B+light_sem_B)+2700,
            (light_mean_B+light_sem_B)+2700, (light_mean_B+light_sem_B)+2200]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1
    )
    
# Light bar 2 to light bar 3 p-value bracket
x_coords = [0.574, 0.574, 0.9, 0.9]
y_coords = [(light_mean_B+light_sem_B)+2500, (light_mean_B+light_sem_B)+2700,
            (light_mean_B+light_sem_B)+2700, (light_mean_C+light_sem_C)+100]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1
    )
    
# Light bar 1 to light bar 3 p-value bracket
x_coords = [0.233, 0.233, 0.9, 0.9]
y_coords = [(light_mean_A+light_sem_A)+4100, (light_mean_A+light_sem_A)+4500,
            (light_mean_A+light_sem_A)+4500, (light_mean_C+light_sem_C)+2900]
for i in range(1,len(x_coords)):
    fig.add_shape(
        type="line",
        xref="paper",
        x0=x_coords[i-1], 
        y0=y_coords[i-1], 
        x1=x_coords[i], 
        y1=y_coords[i],
        line=dict(color='rgba(0,0,0,1)', width=1.5), opacity=1        
    )

    
# P-values
# Darl to Dark
# p-value: Dark bar 1 to dark bar 2
fig.add_annotation(text="p = 0.5111",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.23, y=0.640, showarrow=False,
                   font=dict(size=12, color="black")
                  )

# p-value: Dark bar 1 to dark bar 3
fig.add_annotation(text="p = 0.5898",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.48, y=0.845, showarrow=False,
                   font=dict(size=12, color="black")
                  )

# p-value: Dark bar 2 to dark bar 3
fig.add_annotation(text="p = 0.9910",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.62, y=0.640, showarrow=False,
                   font=dict(size=12, color="black")
                  )
# Light to Light
# p-value: Light bar 1 to light bar 2
fig.add_annotation(text="p = 0.4897",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.325, y=0.75, showarrow=False,
                   font=dict(size=12, color="black")
                  )

# p-value: Light bar 1 to light bar 3
fig.add_annotation(text="p = 0.5841",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.56, y=0.95, showarrow=False,
                   font=dict(size=12, color="black")
                  )

# p-value: Light bar 2 to light bar 3
fig.add_annotation(text="p = 0.9870",
                   name="p-value",                                  
                   xref="paper", yref="paper",
                   x=0.753, y=0.75, showarrow=False,
                   font=dict(size=12, color="black")
                  )




# Customization of layout and traces
fig.update_layout(template='simple_white', title='', height=600, width=1000, yaxis_title='Distance moved [mm]', barmode='group',
                  dragmode='drawrect', font_size=12, hoverlabel_namelength=-1, legend_title_text='',
                  bargroupgap=0,
                 )

fig.update_traces(hoverinfo="x+y")

# Set custom x-axis labels
fig.update_xaxes(title='', tickvals=["A", "B", "C"],
                 ticktext=["Elephant", "Dolphin", "Cat"],
                 ticks="", tickfont_size=14
                )

fig.update_yaxes(range=[0, 7500]
                )


print(df.groupby(['genotype', 'period']).mean())
print(df.groupby(['genotype', 'period']).sem())

fig.show()
   
   
   
   
   
   
   
   
   