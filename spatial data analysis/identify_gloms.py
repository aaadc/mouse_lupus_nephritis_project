import pandas as pd
import numpy as np
import sys
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

import os
os.chdir('')
glom_dat = pd.read_csv('M09909_gloms_position.csv',index_col=0)
coords = np.array(glom_dat[['coor_x','coor_y']])
n_clusters = 95 # different among samples
n_clusters = int(n_clusters)
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(coords)
glom_dat.insert(glom_dat.shape[1], 'cluster_label', kmeans.labels_)


output_path = '/work_path/'
output_prefix = 'sample_id&glom_number'
cluster_label = ''.join([ output_path,output_prefix, "_cluster_label.txt"])
macro_inside_label = ''.join([ output_path,output_prefix, "_macro_inside_label.txt"])
cluster_center = ''.join([ output_path,output_prefix, "_cluster_center.txt"])
cluster_fig = ''.join([ output_path,output_prefix, "_cluster_label.pdf"])
np.savetxt(cluster_label, glom_dat)
np.savetxt(cluster_center, kmeans.cluster_centers_)


plt.scatter(glom_dat['coor_x'],glom_dat['coor_y'],s=3,c=glom_dat['cluster_label'],cmap= 'tab20')
plt.savefig(cluster_fig)
