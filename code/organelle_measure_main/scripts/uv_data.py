# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 11:43:13 2024

@author: Anang
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def normalize_data(bfa):
    
    normalized_avg = {}
    std = {}
    
    for key, value in bfa.items():
        if isinstance(value[0], list):  # Check if the value is a list of lists
            array = np.array(value)
            normalized_array = array / array[:, [0]]  # Normalize each column by its first value
            avg_array = np.mean(normalized_array, axis=0)  # Calculate the average across rows
            
            avg_std = np.std(normalized_array, axis=0)  # Calculate the average across rows
            normalized_avg[key] = avg_array
            std[key]=avg_std
        else:  # If the value is a list of numbers
            normalized_avg[key] = list(map(lambda x: x/value[0],value))
            std[key]=0
    return normalized_avg,std
    
    
    
    
    


control={
    'control':[ [0.199,0.333,0.565,0.901,],
               
               [0.111,0.169,0.297,0.488,],
               
               [0.115,0.182,0.305,0.494,] , 
               [0.073,0.132,0.221,0.323]
        
        ]
    
    }

cntrol_array = np.array(control['control']).T

uv_cnt = cntrol_array / cntrol_array[0]



###chx ug/ml
uv_chx={
        
        'chx:0.01ug/ml ':[],
        
        'chx:0.1ug/ml':[],
        
        'chx:0.5ug/ml':[],
        
        'chx:2ug/ml':[]
        

        
        }

###bfa ug/ml

bfa={
        
        'bfa:25mg/ml':[ [0.116,0.176,0.286,0.474,],
                 
                 [0.186,0.263,0.413,0.668,],
                 
                 [0.105,0.154,0.247,0.413,],
   
            ],
        
 
        'bfa:75mg/ml':[  [0.185,0.227,0.31,0.423,],
                 
                 [0.126,0.174,0.206,0.229,],
                 
                 [0.099,0.123,0.153,0.201,],
            
            ],
        
        'bfa:100mg/ml':[0.185,0.206,0.233,0.249,],

        
        }


uv_bfa,std=normalize_data(bfa)
import matplotlib.pyplot as plt
for key, value in uv_bfa.items():
    plt.plot(value, marker='o', label=key)


plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Uv measurements')
plt.legend()
plt.show()



uv_chxrapa={
    
    'chxrapa: 2:1ug/ml':[1,1.36559,1.45161,1.55914,],
    
    'chxrapa:2:2ug/ml':[1,1.45968,1.68548,1.81452,],
    
    
    }


##rapa ug/ml
uv_rapa={
    
    'rapa:1ug/ml':[1,1.4,1.7,1.87273,],
    
    'rapa:2ug/ml':[1,1.56643,1.95105,2.32867,],
    
    
    }


uv_cu={ 
       
       'control':[ [0.199,0.333,0.565,0.901,],
                  
                  [0.111,0.169,0.297,0.488,],
                  
                  [0.115,0.182,0.305,0.494,] , 
                   
           
           ],
       
       
       
       'cu:2.5 ug/ml':[
           [0.096,0.149,0.205,0.208],
           [0.124,0.202,0.253,0.250],
          
           ],

       'cu:5ug/ml':[
           [0.118,0.207,0.242,0.240],
           [0.104,0.171,0.204,0.198],
           ],
       
       'cu:0.1ug/ml':[
           [0.106,0.134,0.254,0.485]
           ],
       
       'cu:1ug/ml':[
           [0.162,0.345,0.291,0.376]
           
           ],
       
        'cu:0.05ug/ml':[0.110,0.198,0.261,0.284],
       
    

       }




uv_bfa,std=normalize_data(uv_cu)
import matplotlib.pyplot as plt
for key, value in uv_bfa.items():
    plt.plot(value, marker='o', label=key)
    
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Uv measurements')
plt.legend()
plt.show()


uv={
    
    'control':[ [0.199,0.333,0.565,0.901,],
               
               [0.111,0.169,0.297,0.488,],
               
               [0.115,0.182,0.305,0.494,] , 
        ],
    
    
    
    'bfa:25mg/ml':[ [0.116,0.176,0.286,0.474,],
             
             [0.186,0.263,0.413,0.668,],
             
             [0.105,0.154,0.247,0.413,],
        ],
    

    'bfa:75mg/ml':[  [0.185,0.227,0.31,0.423,],
             
             [0.126,0.174,0.206,0.229,],
             
             [0.099,0.123,0.153,0.201,],
        ],
    
    'bfa:100mg/ml':[0.185,0.206,0.233,0.249,],
    
     
    'cu:2.5 ug/ml':[
        [0.096,0.149,0.205,0.208],
        [0.124,0.202,0.253,0.250],
       
        ],

    'cu:5ug/ml':[
        [0.118,0.207,0.242,0.240],
        [0.104,0.171,0.204,0.198],
        ],
    
    'cu:0.1ug/ml':[
        [0.106,0.134,0.254,0.485]
        ],
    
    'cu:1ug/ml':[
        [0.162,0.345,0.291,0.376]
        
        ],

    
    }






##########
marker_styles = {
    'con': 'o',
    'bfa': '^',
    'cu:': 's'
}


# key=uv.keys()




uv_all,std=normalize_data(uv)
# first_few_letters = [key[:3] for key in uv_all]
time=[0,1,2,3]

for key, value in uv_all.items():
    first_few_letters1=key[:3]
    marker = marker_styles.get(first_few_letters1, 'o')  # Get marker style from dictionary, default to 'o' 
    # plt.plot(value, marker=marker, label=key)
    
    plt.errorbar(time,value,yerr = std[f'{key}'], marker=marker, label=key)
    
    
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Uv measurements')
plt.legend()
plt.show()


