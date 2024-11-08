# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 12:09:16 2023

@author: Anang
"""

import scipy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import container
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import organelle_measure.vars_allround1data as var
# from organelle_measure.data import read_results




path_in=var.path_out

folders=[
        "July22",
        "July25"
        ]
subfolder='measurements'


for folder in folders:
    for path_i in Path(f"{path_in}/{folder}/{subfolder}").glob("vacuole*.csv"):
        
        print(path_i)
    

