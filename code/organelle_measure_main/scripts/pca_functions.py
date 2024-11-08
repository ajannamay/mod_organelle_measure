# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 15:42:18 2024

@author: Anang
"""

from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


penguins_data=penguins_data

columns_of_interest = ['Species', "Culmen Length (mm)", "Culmen Length (mm)", 
                       "Flipper Length (mm)", "Body Mass (g)", "Sex"]
penguins_df = penguins_raw.loc[:,columns_of_interest]


# shorten penguin specie name
penguins_df[['Species']]=penguins_df.Species.str.split(" ",expand=True,).loc[:,0]
# replace "." to missing value
penguins_df=penguins_df.replace(".", np.nan)
# drop all rows containing missing value
penguins_df=penguins_df.dropna()

penguins_data=penguins_df.select_dtypes(np.number)
penguins_info=penguins_df.select_dtypes(exclude='float')
penguins_info.Species.unique()

sex=penguins_info.Sex.tolist()
species=penguins_info.Species.tolist()


pca = PCA(n_components=4)
penguins_pca= pca.fit_transform(penguins_data)


pc_df = pd.DataFrame(data = penguins_pca , 
        columns = ['PC1', 'PC2','PC3', 'PC4'])

pc_df['Sex']=sex
pc_df['Species']=species

pca.explained_variance_ratio_


import seaborn as sns
plt.figure(figsize=(12,10))
with sns.plotting_context("notebook",font_scale=1.25):
    sns.scatterplot(x="PC1", y="PC2",
                    data=pc_df, 
                    hue="Species",
                    style="Sex",
                    s=100)
    
###PCA with scaled data

random_state = 0
pca_scaled = make_pipeline(StandardScaler(),
                    PCA(n_components=4, random_state=random_state))

penguins_pc_scaled=pca_scaled.fit_transform(penguins_data)

pca_scaled.named_steps['standardscaler'].fit_transform(penguins_data)




pc_scaled_df = pd.DataFrame(data = penguins_pc_scaled , 
        columns = ['PC1', 'PC2','PC3', 'PC4'])

pca_scaled.named_steps['pca'].explained_variance_ratio_*100


pc_scaled_df['Species'] = species
pc_scaled_df['Sex'] = sex
pc_scaled_df.head()


plt.figure(figsize=(12,10))
with sns.plotting_context("talk",font_scale=1.25):
    sns.scatterplot(x="PC1", y="PC2",
                    data=pc_scaled_df, 
                    hue="Species",
                    style="Sex",
                    s=100)
    plt.xlabel("PC1: "+f'{var_explained[0]:.0f}'+"%")
    plt.ylabel("PC2: "+f'{var_explained[1]:.0f}'+"%")
    plt.title("PCA with Scaled Penguins Data")
plt.savefig("PCA_plot_PC1_vs_PC2_Penguins_scaled_data.png",
                    format='png',dpi=150)


##how much of the variation captured by PC1 is due to Species

sns.plotting_context("talk",font_scale=1.25)
plt.figure(figsize=(8,6))
sns.boxplot(x="Species",y="PC1",
            width=0.5,
            data=pc_scaled_df)
plt.xlabel("Species", size=14)
plt.ylabel("PC1: "+f'{var_explained[0]:.0f}'+"%", size=14)
plt.title("PC1 vs Species: Penguins Data", size=16)
plt.savefig("PCA_plot_PC1_vs_Species_Penguins_scaled_data.png",
                    format='png',dpi=150)


'''it is an indication that all the original variables have some weight in the principal components.

MOre  the angle less the corelation between the variables.


'''
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
biplot(components, np.transpose(PCA.components_), list(iris.columns))