import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

import nglview as nv
import warnings

from options import Options
from plot import Plot

#https://userguide.mdanalysis.org/stable/examples/analysis/reduced_dimensions/pca.html
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')
#%matplotlib inline

try:
    #create universe
    u = mda.Universe("md.gro","md_noPBC.xtc")

    print(mda.Universe("md.gro","md_noPBC.xtc"))
except FileNotFoundError:
    print('Could not find files with extension .pdb .xtc or .gro, please verify inputs.')
    sys.exit()

#align to specific kind of structure
aligner = align.AlignTraj(u, u, select='backbone',
                          in_memory=True).run()

#define pc as PCA to backbone
pc = pca.PCA(u, select='backbone',
             align=False, mean=None,
             n_components=None).run()

opt = Options(pc, u)
plot = Plot()

number_backbone, backbone = opt.pca_analysis()

cumulated_variance = opt.calculate_cumulated_variance(0, 10)
print(f'Cumulated variance: {cumulated_variance}')
plot.cumulated_variance_plot(cumulated_variance, pc, 0, 10)

variance = opt.calculate_variance(0, 10)
print(f'Variance: {variance}')
plot.variance_plot(variance, pc, 0, 10)

combine_df, transformed = opt.combine_dataframe(backbone, 2)
combine_plot = plot.time_plot(combine_df)
combine_plot_scatter = plot.time_scatterplot(combine_df)

visualize = opt.visualize_movement(transformed, backbone)