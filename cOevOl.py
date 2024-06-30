#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 21:27:21 2023

@author: elbangue
"""
import concurrent.futures
#import dash
#import dash_core_components as dcc
#import dash_html_components as html
#import plotly.express as px
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import manifold
import os
import cosa as kosa
import time
#import ray
import psutil
#import sys
import tkinter as tk
#import logging


NORM_FONT = ("Helvetica", 10)

def popupmsg(msg):
    
    popup = tk.Tk()
    popup.wm_title(" PlanTiNT_0.0.ß_emilia ")
    label = tk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = tk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    #popup.mainloop()
    

popupmsg('       Prediction of Protein-Protein Interaction by Coevolution \n Plant_Biology| Saarland University | alexander.banguelacastillo@uni-saarland.de|\n'  )

                 

#cwd = os.getcwd()  
#Call the first two Arms

#cambio para hongo x plantas
#import after_blast_order_to_PPIFT as ABOT
#os.chdir(os.path.dirname(sys.argv[0]))

import PPIDFT_aux_functions_ABC2020_nuevo as ppi
#
#Third Arm start here: modification of PPIFT analysis from Yin, C. & Yau, Stephen S.-T. (2017)
# 


print('\n                 <======================================================================================')
print('| Third_arm_strike       |  PPIFT analysis |    Please wait!  Based on  Yin, C. & Yau, Stephen S.-T. (2017)                                |')
print('                   ======================================================================================>')



# Existing paralela function

def paralela(proteinsA_like, p, files_prot, proteinNames):
    nameB = proteinNames[p]
    proteinsB = ppi.getAllSequences(nameB, files_prot)
    dist_like = 1 - ppi.scorePPITreesP(proteinsA_like, proteinsB)
    return dist_like

# Rest of the code

provisonal = r'/home/elbangue/Desktop/Saarland/Query/Results/Results_BRUT/Results_V'
proteinNames = kosa.files_to_PPIFDT(provisonal)
n = len(proteinNames)
distM = np.zeros([n, n], dtype=np.float128)
distV = []
markerSizes = []
df = pd.DataFrame([])
num_cpus = psutil.cpu_count(logical=False)

with concurrent.futures.ThreadPoolExecutor() as executor:
    for i in range(n):
        start = time.time()
        nameA = proteinNames[i]
        proteinsA = ppi.getAllSequences(nameA, provisonal)
        markerSizes.append(len(proteinsA[0]))

        results = list(executor.map(paralela, [proteinsA]*n, range(n), [provisonal]*n, [proteinNames]*n))

        c = 0
        for x in results:
            distV.append(x)
            distM[i, c] = x
            c += 1

        time.sleep(0.001)
        end = time.time()

        print("\r Loading ... {}".format(i + 1), " out of ", n, " in ", round(end - start, 1), "sec ->", end=nameA)
        data = pd.Series(distV, index=proteinNames)
        distV = []
        df = df.append(data, ignore_index=True)
    
df.index = proteinNames
df = df.sort_index()
df.to_csv(r'/home/elbangue/Desktop/Saarland/Query/Results/dist_file.csv', sep ='\t')
'''
Multidimensional Scaling (MDS) analysis of disstances of PPIs. It maps distance matrix to 2-D
space for the PPI visualization.
'''
# Dimension is 2, the MDS fitting may be different from each run.
mds=manifold.MDS(n_components=2, max_iter=500,dissimilarity="precomputed", n_jobs=1)
pos=mds.fit(distM).embedding_
#planta x hongo
#os.mkdir (ABOT.internal_path +r'/PlanTInt_Final_Results')
os.mkdir (r'/home/elbangue/Desktop/Saarland/Query/Results/PlanTInt_Final_Results')

# matrix write to csv file

#plant x hongo
#f = open(ABOT.internal_path +r'/PlanTInt_Final_Results/PlanTInt_query_proteins_table.csv', 'w')
f = open(r'/home/elbangue/Desktop/Saarland/Query/Results/PlanTInt_Final_Results/PlanTInt_query_proteins_table.csv', 'w')
i =0
for item in pos:
    f.write(proteinNames[i]+','+','.join([str(x) for x in item]) + '\n')
    i = i+1
f.close()
#labels=proteinNames

# Plot the points
# markerSizes= [300]#[676,739,251,289,341,326,1210] #Illustration for realative lengths of amino acid sequences
v=[0.02*j for j in markerSizes]
#colors= ['r'] #['r','b','k','c','y','m','g'] 

for i in range(len(pos)):
    colors = np.random.rand(3,)
    plt.plot(pos[i, 0], pos[i, 1],marker='o',c=colors,markersize=v[i])
    plt.text(pos[i, 0]+0.008, pos[i, 1]+0.008, proteinNames[i], fontsize=3) # a little off the point positions

plt.xlabel('Relative interaction', labelpad=5,fontsize=11)
plt.ylabel('Relative interaction',labelpad=5,fontsize=11)
#plata x hongo
#plt.savefig(ABOT.internal_path +r'/PlanTInt_Final_Results/PlanTInt_query_proteins_figure.png', dpi=1000)
plt.savefig(r'/home/elbangue/Desktop/Saarland/Query/Results/PlanTInt_Final_Results/PlanTInt_query_proteins_figure.png', dpi=1000) 
plt.show()
print ('--------------------------->')
print ('Third arm strike completed!!')
#plata x hongo
#print ('PlanTInt has finished. Look for results in:', ABOT.internal_path,'/PlanTInt_Final_Results')
print ('PlanTInT_0.0.ß_emilia has finished. Look for results in: /home/elbangue/Desktop/Saarland/Query/Results/PlanTInt_Final_Results')
os.system ('ntfy -b telegram send final_procesamiento_PlantInt')