# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:20:41 2024

@author: Anang
"""



##################extra##########################


###=====================================================
  #%%[3] start reading files       
df_bycell_chx=pd.read_csv("../../../Analysis/Analysis_files/df_bycell_Oct10_2023.csv")

an.get_details(exp_names, df_bycell_chx)

df_bycell_old = df_bycell_chx[~df_bycell_chx['folder'].isin(subfolders_chx)]

an.get_details(exp_names, df_bycell_old)

##======================read after re-analysis new files============================
df_bycell_new = read_results(Path(out_chx),subfolders_chx,(px_x,px_y,px_z))

df_bycell_new['experiment'] = df_bycell_new .apply(an.custom_function, axis=1)

df_bycell_concated=pd.concat([df_bycell_old,df_bycell_new ]).reset_index(drop=True)

an.get_details(exp_names,df_bycell_concated)

data=df_bycell_concated.drop(["Unnamed: 0"],axis=1)

#%%data having cell area <700. reduced the data 
data=data[data["cell-volume"]<700]

df_bycell=data



###=== New Data frame======

data1[:,"ratio"]=data(data["mean","total","total-fration"])/data(data["total-fraction"] if data["organelle"]=="ER")


data=data[data["experiment"]!="cycloheximide"]

data["experiment"].unique()

###percentage change calculation
data1 = data.groupby(['condition','organelle']).mean().reset_index()

org=org_names[4]
cond=0.0
prop="total-fraction"
cond1=100.
x0=data1[(data1["organelle"]==f"{org}") & (data1["condition"]==cond)][f"{prop}"].values
x1=data1[(data1["organelle"]==f"{org}") & (data1["condition"]==cond1)][f"{prop}"].values

print(f"percentage chnage:{org}_{prop}",(x0-x1)*100/x0)



#########
result_df = pd.DataFrame(columns=['organelle', 'mean'])

# Create a FacetGrid
g = sns.FacetGrid(data=data, col="organelle",
                  hue="organelle",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False, margin_titles=True)

# Define a function to extract mean values and append to result_df
def extract_means(*args, **kwargs):
    organelle = kwargs['hue']
    mean_value = kwargs['y'].mean()
    result_df.loc[len(result_df)] = [organelle, mean_value]

# Map the pointplot and extract means
g.map(
    extract_means, "growth", "mean",
    dodge=.4, linestyle="none", errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)






#######===mean
data["experiment"].unique()


data1=data[data["growth"]!=0.16]


data1=data1.sort_values(by=['growth','experiment'])

data1['growth'].unique()

xx=data1[data1["experiment"]=="control"]



g = sns.FacetGrid(data=data, col="organelle",
                  # hue="organelle",
                  hue="experiment",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True,
                    despine=False
                  )
g.map(
    sns.pointplot, "growth","mean",
      linestyle="-",errorbar=("ci"),
     
    markersize=20, markeredgewidth=3,
    err_kws={'elinewidth': 1} 
    
)
g.add_legend()





data.info()

g = sns.FacetGrid(data=data1, col="organelle",
                  hue="experiment",
                  # hue_order=experiment_order,
                  # palette=custom_palette,
                  col_wrap=3,
                    sharex=False, 
                  sharey=False,
                  margin_titles=True)

g.map(
    sns.pointplot, "growth", "mean",
    linestyle="-", errorbar=("ci"),
    marker="-",#["^", "o", "*"],
    markersize=20, markeredgewidth=3,
    native_scale=False,
    ax=True,
    err_kws={'elinewidth': 1} 
)

# g.set(xlim=(0, 1))

# Add legend
g.add_legend(title="Experiment")




###########====mean

g = sns.FacetGrid(data=data, col="organelle",
                  hue="experiment",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)

g.map(
    sns.pointplot, "growth","mean",
   
      linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
   
)

###########====count
g = sns.FacetGrid(data=data, col="organelle",
                  hue="experiment",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(
    sns.pointplot, "growth","count",
   
    dodge=.4, linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)

###########====total
g = sns.FacetGrid(data=data, col="organelle",
                  hue="organelle",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(
    sns.pointplot, "growth","total",
   
    dodge=.4, linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)
########== total-fraction
g = sns.FacetGrid(data=data, col="organelle",
                  hue="organelle",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(
    sns.pointplot, "growth","total-fraction",
   
    dodge=.4, linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)



import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from organelle_measure.data import read_results
import organelle_measure.vars_allround1data as var
from termcolor import colored
import organelle_measure.an_functions as an

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from sklearn.decomposition import PCA


def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)


rng = np.random.RandomState(1)
X = np.dot(rng.rand(2, 2), rng.randn(2, 200)).T
pca = PCA(n_components=2, whiten=True)
pca.fit(X)

fig, ax = plt.subplots(1, 2, figsize=(16, 6))
fig.subplots_adjust(left=0.0625, right=0.95, wspace=0.1)

# plot data
ax[0].scatter(X[:, 0], X[:, 1], alpha=0.2)
for length, vector in zip(pca.explained_variance_, pca.components_):
    print(vector)
    
    
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca.mean_, pca.mean_ + v, ax=ax[0])
ax[0].axis('equal');
ax[0].set(xlabel='x', ylabel='y', title='input')

# plot principal components
X_pca = pca.transform(X)
ax[1].scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.2)
draw_vector([0, 0], [0, 3], ax=ax[1])
draw_vector([0, 0], [3, 0], ax=ax[1])
ax[1].axis('equal')
ax[1].set(xlabel='component 1', ylabel='component 2',
          title='principal components',
          xlim=(-5, 5), ylim=(-3, 3.1))

# fig.savefig('figures/05.09-PCA-rotation.png')







def plot_pca_components(x, coefficients=None, mean=0, components=None,
                        imshape=(8, 8), n_components=8, fontsize=12,
                        show_mean=True):
    if coefficients is None:
        coefficients = x
        
    if components is None:
        components = np.eye(len(coefficients), len(x))
        
    mean = np.zeros_like(x) + mean
        

    fig = plt.figure(figsize=(1.2 * (5 + n_components), 1.2 * 2))
    g = plt.GridSpec(2, 4 + bool(show_mean) + n_components, hspace=0.3)

    def show(i, j, x, title=None):
        ax = fig.add_subplot(g[i, j], xticks=[], yticks=[])
        ax.imshow(x.reshape(imshape), interpolation='nearest')
        if title:
            ax.set_title(title, fontsize=fontsize)

    show(slice(2), slice(2), x, "True")
    
    approx = mean.copy()
    
    counter = 2
    if show_mean:
        show(0, 2, np.zeros_like(x) + mean, r'$\mu$')
        show(1, 2, approx, r'$1 \cdot \mu$')
        counter += 1

    for i in range(n_components):
        approx = approx + coefficients[i] * components[i]
        show(0, i + counter, components[i], r'$c_{0}$'.format(i + 1))
        show(1, i + counter, approx,
             r"${0:.2f} \cdot c_{1}$".format(coefficients[i], i + 1))
        if show_mean or i > 0:
            plt.gca().text(0, 1.05, '$+$', ha='right', va='bottom',
                           transform=plt.gca().transAxes, fontsize=fontsize)

    show(slice(2), slice(-2, None), approx, "Approx")
    return fig


from sklearn.datasets import load_digits

digits = load_digits()
sns.set_style('white')

fig = plot_pca_components(digits.data[10],
                          show_mean=False)

pca = PCA(n_components=8)
Xproj = pca.fit_transform(digits.data)
sns.set_style('white')
fig = plot_pca_components(digits.data[10], Xproj[10],
                          pca.mean_, pca.components_)

###choosing number of PC
pca = PCA().fit(digits.data)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance');


data=data
# Compute the mean
m = np.array([data['mean'].mean(), 
              data['count'].mean(),data['total'].mean(),data['cell-volume'].mean(),
              data['total-fraction'].mean()],
             )

# Plot petal length vs petal width only

plt_x='cell-volume'
plt_y='total-fraction'


for experiment, group in data.groupby(['experiment']):
    plt.plot(group[f'{plt_x}'], group[f'{plt_y}'],
               label=experiment, marker='o', linestyle='none')

# Add the mean value to the plot
plt.plot(m[0], m[1],
         marker='*', color='black', markersize=20,
         linestyle='none', label='mean')

plt.legend(loc=0, fontsize=15)
plt.margins(0.01)
plt.xlabel(f'{plt_x}')
plt.ylabel(f'{plt_y}');


###PCA
# Substract the mean from the measurements.
df_cent = data.loc[:, ['mean', 'count','total','cell-volume','total-fraction']]
for col in df_cent.columns:
    df_cent[col] -= df_cent[col].mean()

# Take a look
df_cent.head()
cov_mat = np.cov(m=df_cent.transpose())
print('Covariance matrix \n%s'%cov_mat)

eig_vals, eig_vecs = np.linalg.eig(cov_mat)

print('Eigenvectors\n', eig_vecs)
print('\nEigenvalues\n', eig_vals)


# Plot Petal length vs petal width only
for key, group in data.groupby(['experiment']):
    plt.plot(group[f'{plt_x}'], group[f'{plt_y}'],
               label=key, marker='o', linestyle='none')

# Add the mean value to the plot
plt.plot(m[0], m[1], marker='*', color='black', markersize=20)

# Add arrows showing the eigenvectors
plt.quiver([m[0]]*2, [m[1]]*2, eig_vecs[:,1], eig_vecs[:,0], zorder=11, 
           width=0.01, scale=6)
    
# Tidy up plot
plt.legend(loc=0, fontsize=15)
plt.margins(0.01)
plt.xlabel(f'{plt_x}')
plt.ylabel(f'{plt_y}')


# Compute how much variance is explained by each principal component
print("""
PCA 1: {0:.2f}% of the variance
PCA 2:  {1:.2f}% of the variance
""".format(*tuple(eig_vals / eig_vals.sum() * 100)))

# Project data to our 1D space
df_1D = pd.DataFrame(data[['mean', 'count','total','cell-volume','total-fraction']].dot(eig_vecs[:,0]),
                    columns=['projection'])

# Add back the species column
df_1D['experiment'] = data['experiment']
df_1D.head()


for key, group in df_1D.groupby(['experiment']):
    
    plt.plot(group['projection'], np.ones(len(group)) * 3, alpha=0.7,
             label=key, marker='o', linestyle='none')
    group_num = (np.where(data['experiment'] == key)[0])
    plt.plot(group['projection'], np.ones(len(group)) * group_num, alpha=0.7,
             marker='o', linestyle='none')

plt.margins(0.05)
plt.yticks(range(4), np.append(data['experiment'], 'all'))
plt.xlabel('PCA 1')
#============================================

from sklearn.decomposition import PCA
features = ['mean','count','total','cell-volume','total-fraction']

X = data[features]
# normalizing features
X_norm = (X - X.mean(axis=0))/X.std(axis=0)

# principal component analysis on features
pca = PCA()

# fit and transform X_norm to PCA dataframe
X_pca = pca.fit_transform(X_norm)


# converting to dataframe
names = [f"PC{i+1}" for i in range(X_pca.shape[1])]
X_pcadf = pd.DataFrame(X_pca, columns=names)

X_pcadf['experiment']=data["experiment"]




g = sns.FacetGrid(data=X_pcadf,hue="experiment")
g.map(sns.scatterplot, "PC1",
"PC2",palette='Set1' )
g.add_legend()


plt.scatter(X_pcadf["PC3"],X_pcadf["PC4"])



print(X_pcadf.head())
print("+++++++++++++++++++++++++++++++++++++++++++++++++++")
print("shape of pca df:", X_pcadf.shape)

pca.singular_values_




from sklearn.decomposition import PCA
non_organelle=False
has_volume=True
is_normalized=True
property1="total-fraction"

org_name=["ER","LD","golgi","mitochondria","peroxisome","vacuole"]
columns1 = [*org_name,"non-organelle"] if non_organelle else org_name
num_pc = 7 if non_organelle else 6

df_pca = data.pivot_table(index=["experiment", "condition", "field", "idx-cell", "folder","cell-volume"],
                          columns="organelle", values="total-fraction", aggfunc="first").reset_index()

if is_normalized:
    for col in columns1:
        df_pca[col] = (df_pca[col]-df_pca[col].mean())/df_pca[col].std()

df_pca.reset_index(inplace=True)

# Get Principal Components (PCs)
np_pca = df_pca[columns1].to_numpy()

pca.fit(np_pca).transform(np_pca)
pca_components = pca.components_
pca_var_ratios = pca.explained_variance_ratio_

# pca = PCA(n_components=num_pc)

from sklearn.preprocessing import StandardScaler
# Scale data before applying PCA
scaling=StandardScaler()
df1=np_pca 
# Use fit and transform method 
scaling.fit(df1)
Scaled_data=scaling.transform(df1)
 
# Set the n_components=3
principal=PCA(n_components=3)
principal.fit(Scaled_data)
x=principal.transform(Scaled_data)

x1 = pd.DataFrame(x, columns =['PC1', 'PC2', 'PC3'])


x2=pd.concat([df_pca,x1 ],axis=1).reset_index(drop=True)
 
# Check the dimensions of data after PCA
print(x.shape)
print(df_pca.shape)
print(x2.shape)

plt.figure(figsize=(10,10))
plt.scatter(x[:,0],x[:,1],c=df_pca['experiment'],cmap='plasma')
plt.xlabel('pc1')
plt.ylabel('pc2')

fig = plt.figure(figsize=(10,10))

# choose projection 3d for creating a 3d graph
axis = fig.add_subplot(111, projection='3d')
axis.scatter(x[:,0],x[:,1],x[:,2], cmap='plasma')
axis.set_xlabel("PC1", fontsize=10)
axis.set_ylabel("PC2", fontsize=10)
axis.set_zlabel("PC3", fontsize=10)



g = sns.FacetGrid(data=x2, col="condition",
                  hue="experiment",
                  )
g.map(
    sns.scatterplot, "PC1","PC2",
   
   
)





#%%Plot from using my functions
an.cell_volume_notch(extremes,data)

var="mean"
exp_names=["cycloheximide","control"]
an.get_independent_measurement(org_names,exp_names,df_bycell,var,with_error_bar=False)


var="total-fraction"
exp_names=["cycloheximide","control"]
an.get_independent_measurement(org_names,exp_names,df_bycell,var,with_error_bar=True)


var="total"
exp_names=["cycloheximide","control"]
an.get_independent_measurement(org_names,exp_names,df_bycell,var,with_error_bar=True)


var="count"
exp_names=["cycloheximide","control"]
an.get_independent_measurement(org_names,exp_names,df_bycell,var,with_error_bar=True)


var="mean"
an.get_count_measurement(organelles,exp_names,df_bycell,var,with_error_bar=False)











######========================================================
avg=data.groupby(by=["organelle","condition"]).mean().reset_index()
g = sns.FacetGrid(data=avg, col="organelle",col_wrap=3,sharex=False, sharey=False,margin_titles=True)
g.map(sns.scatterplot, "growth",
"total-fraction",palette='Set1' )
g.add_legend()


data1=data[data["organelle"]=="ER"]

####==============density plots========================

sns.kdeplot(data['cell-volume'], fill=True)

#%%cells and conditions

# data1=data.loc[data['experiment']!="bfa"]

sns.kdeplot(data=data,x="total-fraction",col="organelle",hue="condition",palette="Set1",common_norm=False, bw_method=0.2)

sns.histplot(data=data,x="cell-volume",hue="condition",palette="Set1",kde=True,fill=False)

sns.displot(data=data,x="total-fraction",hue="condition",kde=True,palette="Set1",multiple="stack", fill=False,col="organelle",col_wrap=3)

plt.xlim(-100, 1000)


#%%plots using sns library

g = sns.FacetGrid(data=data, col="organelle", 
                  hue="condition",col_wrap=3,
                  sharex=False, sharey=False,margin_titles=True)
g.map(sns.scatterplot, "cell-volume",
"total-fraction",palette='Set1' )
g.add_legend()

####=====Cell properties=======

##

g = sns.FacetGrid(data=data, col="experiment", 
                  
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(
    sns.pointplot, "condition","cell-volume",
   
    dodge=.4, linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)

###===================================================
g = sns.FacetGrid(data=data, col="organelle", 
                  hue="organelle",col_wrap=3,
                  sharex=False, sharey=False,margin_titles=True)
g.map(sns.scatterplot, "cell-volume",
"count",palette='Set1' )
g.add_legend()
##########################

####=====Ratio plot ========


##################################################
avg=data.groupby(by=["count","organelle","condition","folder"]).mean().reset_index()
g = sns.FacetGrid(data=avg, col="organelle", 
                  hue="condition",
                  col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(sns.scatterplot, "cell-volume",
"total-fraction",palette='Set1' )
g.map(plt.errorbar, "C", "D",yerr=0.5, fmt='o')
g.add_legend()



sns.pointplot(
    data=data, x="condition", y="total-fraction", hue="organelle",
    col="organelle",
    dodge=.4, linestyle="none", errorbar=("pi", 100),
    marker="_", markersize=20, markeredgewidth=3,
)



#####################################################

sns.displot(data=data, x="cell-volume", y="total-fraction")

sns.displot(data=data, x="cell-volume", y="total-fraction",kind="kde", rug=True)

sns.ecdfplot(data=data, x="count")


sns.kdeplot(data=data, x="cell-volume")
sns.rugplot(data=data, x="cell-volume")


sns.scatterplot(data=data, x="cell-volume", y="total-fraction", s=5)
sns.rugplot(data=data, x="cell-volume", y="total-fraction", lw=1, alpha=.005)

sns.displot(data=data,x="cell-volume",kind="kde")

####=============some t-test======================================
from scipy.stats import ttest_ind

dataset1=data[data["condition"].eq(0.0)]["cell-volume"]
dataset2=data[data["condition"].eq(2.)]["cell-volume"]#[:len(dataset1)]

ttest_ind(dataset1, dataset2)

from scipy.stats import mannwhitneyu

mannwhitneyu(dataset1, dataset2)
from scipy.stats import ks_2samp

ks_2samp(dataset1, dataset2)

from scipy.stats import chi2_contingency ###need to be same size data set

data = [dataset1, dataset2]
stat, p, dof, expected = chi2_contingency(data)

import scipy.stats as st 


dd=data[data["condition"].eq(0.01)]
dd["cell-volume"].mean()
dd["cell-volume"].std()

st.t.interval(alpha=0.95, 
              df=len(dd["cell-volume"])-1, 
              loc=np.mean(dd["cell-volume"]),  
              scale=st.sem(dd["cell-volume"])) 


from scipy.stats import norm
import numpy as np
std_dev1 = np.std(dataset1,ddof=1)
std_dev2 = np.std(dataset2,ddof=1)

# Calculate the z-statistic
mean1 = np.mean(dataset1)
mean2 = np.mean(dataset2)
n1 = len(dataset1)
n2 = len(dataset2)
z_statistic = (mean1 - mean2) / np.sqrt((std_dev1**2 / n1) + (std_dev2**2 / n2))

# Calculate the p-value
p_value = 2 * (1 - norm.cdf(np.abs(z_statistic)))

# Print the results
print("Z-statistic:", z_statistic)
print("P-value:", p_value)


from statsmodels.stats.weightstats import ztest as ztest
ztest(dataset1, dataset2, value=0) 



####===========correlation matrix=========================================================

an.corelation_plot(extremes,df_corrcoef,columns,properties,cmap=False)




#####====================================Volume cells===============================================================
#####===================================================================================================



values=sorted(list(df_bycell["condition"].unique()))
   
value_to_number = {value: number for number, value in enumerate(values, 1)}
# sub=["hu","bfa"]
# df=data[~data['experiment'].isin(sub)]

fig, ax = plt.subplots(1, 1)
for condition in data["condition"].unique():
    key=value_to_number[condition]
    # key=condition
    print(value_to_number,key)
    xx=data[data["condition"]==condition]
    yy=xx["cell-volume"].mean()
    yyer=xx["cell-volume"].std()
    plt.errorbar(key,yy,yerr=yyer,fmt="*",ecolor="blue")
    plt.scatter(key,yy,marker="*",s=500,alpha=0.5,cmap="blue") #cmap="nipy_spectral
    

    plt.text(key-0.1,yy,f"{condition}",color="blue",fontsize=12)
     
    if condition ==0:
        xx=data[data["condition"]==condition]
        val=xx["cell-volume"].mean()
        plt.axhline(y = val, color = 'gray', label = 'axvline - full height')
        print(f"experiment: {condition}", val)



#####===================================================================================================

            
org_mean=data.groupby(by=["organelle","condition","count"]).mean().reset_index()
# org_std=df_bycell.groupby(by=["condition","organelle","count","mean"]).std().reset_index()

var_x="cell-volume"
var_y="count"

for org in organelles:
    fig, ax = plt.subplots(1, 1)
    for condition in org_mean["condition"].unique():
        df_bycell1 = org_mean[(org_mean["organelle"] == org) & (org_mean["condition"] == condition)]
        # df_bycell2 = org_mean[(org_std["organelle"] == org) & (org_std["condition"] == condition)]
        
        
        # plt.errorbar(df_bycell1["count"],df_bycell1["mean"],yerr=df_bycell2["mean"])
        
        ax.scatter(df_bycell1[f"{var_x}"],df_bycell1[f"{var_y}"],label=f"Condition: {condition}")
        
        
        ax.legend(loc='upper right',fontsize=font_size)
        
    plt.xlabel(f"{var_x} ",fontsize=font_size)
    plt.ylabel(f"{var_y}_{org}",fontsize=font_size)
    # plt.savefig(f"{chx_im}/sizNnumber_{org}_{condition}.png")
    plt.show()

#####=====================================Scatter plot==============================================================  

            
org_mean=data
# org_std=df_bycell.groupby(by=["condition","organelle","count","mean"]).std().reset_index()

var_x="cell-volume"
var_y="total-fraction"

for org in organelles:
    fig, ax = plt.subplots(1, 1)
    for condition in org_mean["condition"].unique():
        
        df_bycell1 = org_mean[(org_mean["organelle"] == org) & (org_mean["condition"] == condition)]
        # df_bycell2 = org_mean[(org_std["organelle"] == org) & (org_std["condition"] == condition)]
        
        
        # plt.errorbar(df_bycell1["count"],df_bycell1["mean"],yerr=df_bycell2["mean"])
        
        ax.scatter(df_bycell1[f"{var_x}"],df_bycell1[f"{var_y}"],label=f"Condition: {condition}")

        sns.regplot(x=f"{var_x}", y=f"{var_y}", data=df_bycell1)
        plt.title('Regression Plot of X and Y')

        ax.legend(loc='upper right',fontsize=font_size)
        
    plt.xlabel(f"{var_x} ",fontsize=font_size)
    plt.ylabel(f"{var_y}_{org}",fontsize=font_size)
    # plt.savefig(f"{chx_im}/sizNnumber_{org}_{condition}.png")
    plt.show()





#################plots###########

data1=data[["mean","organelle"]]

sns.countplot(x='organelle', data=data1,)
plt.show()

sns.jointplot(x="total-fraction", y="cell-volume", data=data, height=5)
plt.show()



data_chx=data[data["experiment"]!="bfa"]

#################density plot
sns.FacetGrid(data_chx, hue="condition", size=5) \
   .map(sns.kdeplot, "total-fraction") \
   .add_legend()
plt.show() 


df_cor=df_corrcoef.drop(['folder','condition','field','idx-cell','experiment'], axis=1)
sns.heatmap(df_cor.corr(method='pearson'),
            annot = True);
plt.show()

#####============================plots
var="mean"
exp_names=["control"]
an.get_each_date_field(organelles,exp_names,data,var)

exp_names=["control","cycloheximide"]
df_bycell=data



for org in organelles:
    fig, ax = plt.subplots(1, 1)
    data=[]
    df_condition1=[]
    ver_line=-1
    for exp in exp_names:
        print(exp)
        
        
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
        
        
        

###########################
import seaborn as sns
sns.set_theme(style="ticks")

xx_pair=df_bycell_concated.drop(columns=[ 'field','idx-cell','folder','condition'])
sns.pairplot(xx_pair, hue="experiment", height=3, diag_kind="kde")
plt.show() 


df_bycell_concated1=df_bycell_concated[df_bycell_concated["cell-volume"]<5000]

data.plot(kind="scatter", x="total-fraction", y="cell-volume",color="green",s=70 )
plt.show()

data=df_bycell_concated1[["cell-volume","total-fraction","experiment","organelle","condition"]].reset_index(drop=True)

sns.jointplot(x="total-fraction", y="cell-volume", data=data, height=5)
plt.show()

#####=============================================================================================================================
# data=df_bycell_concated1[["cell-volume","total-fraction","organelle"]].reset_index(drop=True)
      


# Modify the graph above by assigning each species an individual color.
sns.FacetGrid(data, hue="condition", size=5) \
   .map(plt.scatter, "total-fraction","cell-volume") \
   .add_legend()
plt.show()       


# To plot the species data using a box plot:

# sns.boxplot(x="total-fraction", y="cell-volume", data=data )
# plt.show()


from pandas.plotting import andrews_curves
andrews_curves(data.drop("organelle", axis=1), "experiment",colormap='rainbow')
plt.show()




####test pca

non_organelle=False
df_bycell_concated1=df_bycell_concated[df_bycell_concated["experiment"]!="bfa"]
columns = [*organelles,"non-organelle"] if non_organelle else organelles
idx = df_bycell_concated1.groupby(["experiment","condition","field","idx-cell","folder"]).count().index
df_pca = pd.DataFrame(index=idx,columns=columns)
for orga in columns:
    print(orga)
    df_pca[orga]=df_bycell_concated1.loc[df_bycell_concated1["organelle"].eq(orga),property]

X= df_pca[columns].to_numpy()

y =  df_pca.experiment

# Standardize the data (mean=0, variance=1)
mean = np.mean(X, axis=0)
std = np.std(X, axis=0)
X = (X - mean) / std












































#%%ifo about the data
data.info()
data.describe()

data.isnull().sum()

data.value_counts("experiment")

#====================================need in analysis for me========================================================================

exp_names=["bfa","control"]

var="mean"
an.an_hist_plot(org_names,exp_names,df_bycell,var)


an.get_each_date_field(organelles,exp_names,data,var)

####==== organelle plot with property

an.plot_cell(data,"cell-volume",growth_rate)




prop="mean"
an.plot_organelle(data,organelles,prop,growth_rate)

an.plot_organelle_field(data,organelles,prop,growth_rate)


####Hist plot cell volume

### sns does bootstrap over 1000 times



sns.pointplot(data=data,x="growth",y="cell-volume")


sns.kdeplot(data=data,x="cell-volume",hue="condition",palette="Set1",common_norm=False, bw_method=0.2)


sns.relplot(data=data,x="growth",y="cell-volume",hue="experiment",style="experiment")

sns.pointplot(data=data,x="growth",y="cell-volume",hue="experiment",style="experiment")
####===========org property======================================================
##################===volume fraction


#%% =====PCA
####====PCA Analysis=====================================###
non_organelle=False
has_volume=True
is_normalized=True
property1="total-fraction"

org_name=["ER","LD","golgi","mitochondria","peroxisome","vacuole"]
columns1 = [*org_name,"non-organelle"] if non_organelle else org_name
num_pc = 7 if non_organelle else 6

####==extract the data based on property here "total-fraction"
df_pca = data.pivot_table(index=["experiment", "condition", "field", "idx-cell", "folder","cell-volume"],
                          columns="organelle", values=f"{property1}", aggfunc="first").reset_index()



df_pca.info()

np_pca = df_pca[columns1].to_numpy()

scaling=StandardScaler()
 
# Use fit and transform method 
scaling.fit(np_pca)
Scaled_data=scaling.transform(np_pca)
 
# Set the n_components=3
principal=PCA(n_components=num_pc)
principal.fit(Scaled_data)
pca_df=principal.transform(Scaled_data)

x1 =pd.DataFrame(pca_df, columns=[f'PC_{i}' for i in range(1, num_pc + 1)])

club_pc_with_data=pd.concat([df_pca,x1 ],axis=1).reset_index(drop=True)


sns.lmplot( x="PC_1", y="PC_2",
  data=club_pc_with_data, 
  fit_reg=False, 
  hue='condition', # color by cluster
  legend=True,
  scatter_kws={"s": 80}) # specify the point size



def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)


rng = np_pca 
X = np_pca
pca = PCA(n_components=6, whiten=True)
pca.fit(X)

fig, ax = plt.subplots(1, 2, figsize=(16, 6))
fig.subplots_adjust(left=0.0625, right=0.95, wspace=0.1)


# plot data
ax[0].scatter(X[:, 0], X[:, 1], alpha=0.2)
for length, vector in zip(pca.explained_variance_, pca.components_):
    print(vector)
    
    
    v = vector * 3 * np.sqrt(length)
    draw_vector(pca.mean_, pca.mean_ + v, ax=ax[0])
ax[0].axis('equal');
ax[0].set(xlabel='x', ylabel='y', title='input')

# plot principal components
X_pca = pca.transform(X)
ax[1].scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.2)
draw_vector([0, 0], [0, 3], ax=ax[1])
draw_vector([0, 0], [3, 0], ax=ax[1])
ax[1].axis('equal')
ax[1].set(xlabel='component 1', ylabel='component 2',
          title='principal components',
          xlim=(-5, 5), ylim=(-3, 3.1))


#####==============PC plot===============
x="PC_1"
y="PC_2"
sns.relplot(
    data=club_pc_with_data, x=f"{x}", y=f"{y}",
    col="experiment", hue="condition", 
    kind="scatter",palette="Set1"
)



pp=sns.pairplot(club_pc_with_data,
             x_vars=["PC_1","PC_2"],
             y_vars=["PC_3","PC_4"],
             hue="condition",palette="Set1")

pp.fig.suptitle(f"PCA using {property1}")


chx_cond=2.
cont_cond=0.
bfa_cond=25.

#####3D PCA plot
PC1_a = club_pc_with_data.loc[club_pc_with_data["condition"] == chx_cond,
                "PC_1"]
PC2_a = club_pc_with_data.loc[club_pc_with_data["condition"] == chx_cond,
                "PC_2"]

PC3_a = club_pc_with_data.loc[club_pc_with_data["condition"] == chx_cond,
                "PC_3"]
 



 
PC1_b =  club_pc_with_data.loc[ club_pc_with_data["condition"] == cont_cond,
                "PC_1"]
PC2_b =  club_pc_with_data.loc[ club_pc_with_data["condition"] ==cont_cond,
                "PC_2"]
PC3_b =  club_pc_with_data.loc[ club_pc_with_data["condition"] ==cont_cond,
                "PC_3"]



PC1_c =  club_pc_with_data.loc[ club_pc_with_data["condition"] == bfa_cond,
                "PC_1"]
PC2_c =  club_pc_with_data.loc[ club_pc_with_data["condition"] ==bfa_cond,
                "PC_2"]
PC3_c =  club_pc_with_data.loc[ club_pc_with_data["condition"] ==bfa_cond,
                "PC_3"]

###2D plot
fig = plt.figure(figsize=(12,7))
ax = fig.add_subplot(111)
 
ax.scatter(PC1_a, 
            PC2_a, 
            c="blue",
           label=f"cycloheximide_{chx_cond}")
 
ax.scatter(PC1_b, 
            PC2_b, 
            c="orange",
           label=f"control_{cont_cond}")

ax.scatter(PC1_c, 
            PC2_c, 
            c="green",
           label=f"bfa_{bfa_cond}")
 
ax.legend(title="Label")
 
plt.title("Figure 1",
          fontsize=16)
plt.xlabel('First Principal Component',
           fontsize=16)
plt.ylabel('Second Principal Component',
           fontsize=16)
#####

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(PC1_a, PC2_a, PC3_a, marker="*")

ax.scatter(PC1_b, PC2_b, PC3_b, marker="^")

ax.scatter(PC1_c, PC2_c, PC3_c, marker="o")
legend1 = ax.legend(['CHX','Control','Bfa'], title="Legend")


ax.set_title("Data Analysis ")
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')



#############################End of PCA Analysis=============
####===t-test

import pandas as pd
import seaborn as sns
from itertools import combinations
from scipy.stats import ttest_ind

# Assume you have four datasets: df1, df2, df3, df4

prop="cell-volume"

df5=data[data["condition"]==0,f"{prop}"]

df1=data.loc[data["condition"]==0,f"{prop}"]
df2=data.loc[data["condition"]==0.01,f"{prop}"]
df3=data.loc[data["condition"]==0.1,f"{prop}"]
df4=data.loc[data["condition"]==2,f"{prop}"]

datasets = [df1, df2, df3, df4]
labels = ["Dataset 1", "Dataset 2", "Dataset 3", "Dataset 4"]

# Create a matrix to store p-values
p_values_matrix = pd.DataFrame(index=labels, columns=labels)

# Perform t-test for all pairs
for pair in combinations(range(len(datasets)), 2):
    idx1, idx2 = pair
    label1, label2 = labels[idx1], labels[idx2]
    t_stat, p_value = ttest_ind(datasets[idx1], datasets[idx2])
    p_values_matrix.at[label1, label2] = p_value*10000
    p_values_matrix.at[label2, label1] = p_value*10000  # The matrix is symmetric
    print(p_value)

# Plot the p-values matrix
plt.figure(figsize=(10, 8))
sns.heatmap(p_values_matrix.astype(float), annot=True, cmap="coolwarm", fmt=".6f")
plt.title('T-Test P-Values Matrix')
plt.show()



####==================corelation matrx plot============================

import plotly.express as px
# Correlation coefficient
for folder in extremes:
    xx=df_corrcoef[df_corrcoef["experiment"]==folder]
    # print(xx["condition"].unique())
    # cond1=str(xx["condition"].iloc[0] )
    np_corrcoef = xx.loc[:,columns].to_numpy()
    corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
    fig = px.imshow(
            corrcoef,
            x=columns,y=columns,
            color_continuous_scale = "RdBu_r",range_color=[-1,1],
            # title=str(xx["condition"].iloc[0] )
        )
    
    fig.write_html(f'{Path(f"{chx_im}/Dec18")}/1corrcoef-nocond_{folder}_all.html')

    for condi in xx["condition"].unique():
        np_corrcoef = xx.loc[df_corrcoef['condition']==condi,['effective-length','cell-area','cell-volume',*properties]].to_numpy()
        corrcoef = np.corrcoef(np_corrcoef,rowvar=False)
        fig = px.imshow(
                corrcoef,
                x=columns,y=columns,
                color_continuous_scale="RdBu_r",range_color=[-1,1]
            )
        fig.write_html(f'{Path(f"{chx_im}/Dec18")}/1corrcoef-nocond_{folder}_{str(condi)}_indivual.html')




####========Extra============
property1="total-fraction"
####==extract the data based on property here "total-fraction"
df_pca = data.pivot_table(index=["experiment", "condition", "field", "idx-cell", "folder","cell-volume"],
                          columns="organelle", values=f"{property1}", aggfunc="first").reset_index()



columns_to_exclude = ['experiment',"condition","field","idx-cell","folder"]

columns_to_divide = [col for col in df_pca.columns if col not in columns_to_exclude]

divisor_column=df_pca["non-organelle"]

# Dividing selected columns by the divisor column
result_df = df_pca[columns_to_divide].div(divisor_column, axis=0)
result_df = pd.concat([df_pca[columns_to_exclude], result_df], axis=1)

result_df=result_df.dropna()

g = sns.FacetGrid(data=result_df, col="experiment",
                   hue="condition",
                  # col_wrap=3,
                  sharex=False, 
                  sharey=False,margin_titles=True)
g.map(
    sns.pointplot, "condition","LD",
   
    dodge=.4, linestyle="none",errorbar=("ci"), 
    marker="_", markersize=20, markeredgewidth=3,
)

sns.pointplot(data=result_df,x="condition",y="peroxisome")



# libraries
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
 
# Get the iris dataset
sns.set_style("white")
df = sns.load_dataset('iris')

# create figure
my_dpi=96
plt.figure(figsize=(480/my_dpi, 480/my_dpi), dpi=my_dpi)
 
# Keep the 'species' column appart + make it numeric for coloring
df['species']=pd.Categorical(df['species'])
my_color=df['species'].cat.codes
df = df.drop('species', 1)
 
# Run The PCA
pca = PCA(n_components=3)
pca.fit(df)
 
# Store results of PCA in a data frame
result=pd.DataFrame(pca.transform(df), columns=['PCA%i' % i for i in range(3)], index=df.index)
 
# Plot initialisation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(result['PCA0'], result['PCA1'], result['PCA2'], c=my_color, cmap="Set2_r", s=60)
 
# make simple, bare axis lines through space:
xAxisLine = ((min(result['PCA0']), max(result['PCA0'])), (0, 0), (0,0))
ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
yAxisLine = ((0, 0), (min(result['PCA1']), max(result['PCA1'])), (0,0))
ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
zAxisLine = ((0, 0), (0,0), (min(result['PCA2']), max(result['PCA2'])))
ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
 
# label the axes
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
ax.set_title("PCA on the iris data set")
ax.legend()
plt.show()

#### corelation matrix 


import seaborn as sns
from matplotlib.pyplot import gcf
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np

# fig, axs = plt.subplots()

networks = sns.load_dataset("brain_networks", index_col=0, header=[0, 1, 2])

# Label 1
network_labels = networks.columns.get_level_values("network")
network_pal = sns.cubehelix_palette(network_labels.unique().size, light=.9, dark=.1, reverse=True, start=1, rot=-2)
network_lut = dict(zip(map(str, network_labels.unique()), network_pal))
network_colors = pd.Series(network_labels, index=networks.columns).map(network_lut)

# Label 2
node_labels = networks.columns.get_level_values("node")
node_pal = sns.cubehelix_palette(node_labels.unique().size)
node_lut = dict(zip(map(str, node_labels.unique()), node_pal))
node_colors = pd.Series(node_labels, index=networks.columns).map(node_lut)

# Label 3
lab3_labels = networks.columns.get_level_values("node")
lab3_pal = sns.color_palette("hls", lab3_labels.unique().size)
lab3_lut = dict(zip(map(str, lab3_labels.unique()), lab3_pal))
lab3_colors = pd.Series(lab3_labels, index=networks.columns, name='lab3').map(lab3_lut)

# Label 4
lab4_labels = networks.columns.get_level_values("node")
lab4_pal = sns.color_palette("husl", lab4_labels.unique().size)
lab4_lut = dict(zip(map(str, lab4_labels.unique()), lab4_pal))
lab4_colors = pd.Series(lab4_labels, index=networks.columns, name='lab4').map(lab4_lut)

network_node_colors = pd.DataFrame(network_colors).join(pd.DataFrame(node_colors)).join(pd.DataFrame(lab3_colors)).join(pd.DataFrame(lab4_colors))

g = sns.clustermap(networks.corr(),
    row_cluster=False, col_cluster=False,
    row_colors = network_node_colors,
    col_colors = network_node_colors,
    linewidths=0,
    xticklabels=False, yticklabels=False,
    center=0, cmap="vlag")


# add legends
for label in network_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=network_lut[label], label=label, linewidth=0);
l1 = g.ax_col_dendrogram.legend(title='Network', loc="center", ncol=5, bbox_to_anchor=(0.35, 0.89), bbox_transform=gcf().transFigure)

for label in node_labels.unique():
    g.ax_row_dendrogram.bar(0, 0, color=node_lut[label], label=label, linewidth=0);
l2 = g.ax_row_dendrogram.legend(title='Node', loc="center", ncol=2, bbox_to_anchor=(0.66, 0.89), bbox_transform=gcf().transFigure)

# create a list for the bar plot patches
xx = []
for label in lab3_labels.unique():
    x = g.ax_row_dendrogram.bar(0, 0, color=lab3_lut[label], label=label, linewidth=0)
    xx.append(x)

# add the legend
legend3 = plt.legend(xx, lab3_labels.unique(), loc="center", title='lab3', bbox_to_anchor=(.78, 0.89), bbox_transform=gcf().transFigure)


# create a list for the bar plot patches
yy = []
for label in lab4_labels.unique():
    y = g.ax_row_dendrogram.bar(0, 0, color=lab4_lut[label], label=label, linewidth=0)
    yy.append(y)

    
# add the second legend
legend4 = plt.legend(yy, lab4_labels.unique(), loc="center", title='lab4', ncol=2, bbox_to_anchor=(.9, 0.89), bbox_transform=gcf().transFigure)
plt.gca().add_artist(legend3)





import seaborn as sns; sns.set(color_codes=True)
import string
iris = sns.load_dataset("iris")
species = iris.pop("species")
lut = dict(zip(species.unique(), "rbg"))
samples = np.repeat(list(string.ascii_letters[0:8]),20)[:150]
sample_cols = dict(zip(set(samples), sns.color_palette("cubehelix", 8)))

row_colors = pd.DataFrame({'species':species.map(lut),
                          'sample':[sample_cols[i] for i in samples]})
g = sns.clustermap(iris, row_colors=row_colors,row_cluster=False)




import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Generate random data
np.random.seed(0)
experiments = ["Exp1", "Exp2", "Exp3"]
growth_types = ["Type1", "Type2", "Type3", "Type4"]

data = pd.DataFrame(np.random.randn(len(experiments), len(growth_types)), columns=growth_types, index=experiments)

# Visualize correlation matrix
correlation_matrix = data.corr()

plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='viridis', fmt=".2f", linewidths=.5)
plt.title('Correlation Matrix of Growth Types')
plt.show()

# Visualize clustered heatmap
sns.clustermap(data, cmap='viridis', method='average', 
               metric='euclidean', col_cluster=True, row_cluster=True)
plt.title('Clustered Heatmap of Experiments and Growth Types')
plt.show()


# Prepare a vector of color mapped to the 'cyl' column

df_corrcoef1=df_corrcoef.loc[:,columns]


my_palette = dict(zip(df_corrcoef.experiment.unique(), ["orange","yellow","brown"]))
row_colors = df_corrcoef.experiment.map(my_palette)
 
# plot
sns.clustermap(df_corrcoef1, metric="correlation", method="single", cmap="Blues", standard_scale=1, row_colors=row_colors)
plt.show()



def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
                      width=ell_radius_x * 2,
                      height=ell_radius_y * 2,
                      facecolor=facecolor,
                      **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


method = PCA(n_components=2, whiten=True)  # project to 2 dimensions
projected = method.fit_transform(np.array(inputs[tags['datum']].tolist()))

figure = pyplot.figure()
axis = figure.add_subplot(111)
# Display data
for label in labels:
  color = np.expand_dims(np.array(settings.get_color(label)), axis=0)
  pyplot.scatter(projected[labels == label, 0], projected[labels == label, 1],
                           c=color, alpha=0.5, label=label, edgecolor='none')

# Centroids
for label in labels:
# Centroids
color = np.array(settings.get_color(label))
# Ellipsis
Views.confidence_ellipse(projected[labels == label, 0], projected[labels == label, 1], axis,
                         edgecolor=color, linewidth=3, zorder=0)


extract_prop = ['experiment','growth','effective-length',*properties]

df_corrcoef1=df_corrcoef.loc[:,extract_prop]


pca_data=df_corrcoef1.select_dtypes(np.number)

pca_data=pca_data.drop(["growth"],axis=1)


labels=pca_data.keys()


NPC=pca_data.shape[1]

boston = pca_data

from sklearn.preprocessing import StandardScaler
x = StandardScaler().fit_transform(boston)
x = pd.DataFrame(x, columns=labels)


from sklearn.decomposition import PCA
pcamodel = PCA(n_components=NPC)
pca = pcamodel.fit_transform(x)
pca.shape

plt.bar(range(1,len(pcamodel.explained_variance_ )+1),pcamodel.explained_variance_ )
plt.ylabel('Explained variance')
plt.xlabel('Components')
plt.plot(range(1,len(pcamodel.explained_variance_ )+1),
         np.cumsum(pcamodel.explained_variance_),
         c='red',
         label="Cumulative Explained Variance")
plt.legend(loc='upper left')

plt.plot(pcamodel.explained_variance_ratio_)
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance')
plt.show()

plt.scatter(pca[:, 2], pca[:, 3])


#Make Plotly figure
import plotly.plotly as py
import plotly.graph_objs as go

# plotly.tools.set_credentials_file(username='prasadostwal', api_key='yourapikey')
# api key hidden

fig1 = go.Scatter3d(x=pca[:, 0],
                    y=pca[:, 1],
                    z=pca[:, 2],
                    marker=dict(opacity=0.9,
                                reversescale=True,
                                colorscale='Blues',
                                size=5),
                    line=dict (width=0.02),
                    mode='markers')

#Make Plot.ly Layout
mylayout = go.Layout(scene=dict(xaxis=dict( title="PCA1"),
                                yaxis=dict( title="PCA2"),
                                zaxis=dict(title="PCA3")),)

#Plot and save html
py.iplot({"data": [fig1],
                     "layout": mylayout},
                     auto_open=True,
                     filename=("3DPlot.html"))




ax = sns.heatmap(pcamodel.components_,
                 cmap='YlGnBu',
                 yticklabels=[ "PCA"+str(x) for x in range(1,pcamodel.n_components_+1)],
                 xticklabels=list(x.columns),
                 annot=True,
                 fmt=".2f",
                 cbar_kws={"orientation": "vertical"})
ax.set_aspect("equal")

def myplot(score,coeff,labels=None):
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

myplot(pca[:,0:2],np.transpose(pcamodel.components_[0:2, :]),list(x.columns))
plt.show()


