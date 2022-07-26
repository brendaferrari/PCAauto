import pandas as pd
import numpy as np
import MDAnalysis as mda
import glob
import csv

class Options:

    def __init__(self, pc, universe):
        self.pc = pc
        self.u = universe

    def pca_analysis(self):

        backbone = self.u.select_atoms('backbone')
        n_bb = len(backbone)
        print(f'There are {n_bb} backbone atoms in the analysis')
        print(f'Shape of pc components: {self.pc.p_components.shape}')

        return n_bb, backbone

    def calculate_variance(self, zeroth_dimension, nth_dimension):

        variance = self.pc.variance[zeroth_dimension:nth_dimension]
        np.savetxt('variance.txt', variance, delimiter=' ')

        return variance

    def calculate_cumulated_variance(self, zeroth_dimension, nth_dimension):

        cumulated_variance = self.pc.cumulated_variance[zeroth_dimension:nth_dimension]
        np.savetxt('cumulated_variance.txt', cumulated_variance, delimiter=' ')

        return cumulated_variance
        
    def calculate_explained_variance(self, pc):

        explained_variance = [('PC1',self.pc.cumulated_variance[0]), ('PC2',pc.cumulated_variance[1]-self.pc.cumulated_variance[0]), 
                              ('PC3',self.pc.cumulated_variance[2]-self.pc .cumulated_variance[1])]

        return explained_variance

    def combine_dataframe(self, backbone, n_dim):

        transformed = self.pc.transform(backbone, n_components=n_dim)
        transformed.shape
        print(f'Transformed shape: {transformed.shape}')

        df = pd.DataFrame(transformed,
                          columns=[f'PC{i+1}' for i in range(n_dim)])
        df['Time (ps)'] = df.index * self.u.trajectory.dt
        df.head()

        df.to_csv('dataFrame.csv', index=True)
        df.to_excel('dataFrame.xlsx', engine='xlsxwriter') 

        return df, transformed

    def visualize_movement(self, transformed, backbone, pc=0):

        print(f'You are calculating tranjectory for {pc+1}')
        pc1 = self.pc.p_components[:, pc]
        trans1 = transformed[:, pc]
        projected = np.outer(trans1, pc1) + self.pc.mean.flatten()
        coordinates = projected.reshape(len(trans1), -1, 3)

        proj1 = mda.Merge(backbone)
        proj1.load_new(coordinates, order="fac")

        proj1.atoms.write('pc'+ str(pc) + '.pdb')
        proj1.atoms.write('pc'+ str(pc) + '.dcd', frames='all')

        return proj1

    def generate_dat_file(self, cumulated_variance=False, variance=True):

        data = []

        import os
        print(os.getcwd())
        #os.chdir('../input')
        workdir = os.getcwd()

        if variance == True:
            File = 'variance.txt'

            for root, dirs, files in os.walk(workdir):
                if File in files:
                    summary_path = root + '/' + File
                    column_name = root.split('/')
                    search_file = open(summary_path, 'r', encoding='utf-8')
                    lines = search_file.readlines()
                    line_clean = []
                    for line in lines:
                        line_cleaned = line.strip('\n')
                        line_clean.append(line_cleaned)
                    line_clean.insert(0, column_name[-1])
                    data.append(line_clean)

            df = pd.DataFrame(data).transpose()
            df.rename(columns=df.iloc[0], inplace=True)
            df = df.iloc[1:]
            df.to_csv('variance.csv', sep=' ', quoting=csv.QUOTE_ALL, index=None)
            df.to_csv('variance.dat', sep=' ', index=None, header=None)

        if cumulated_variance == True:
            File = 'cumulated_variance.txt'

            for root, dirs, files in os.walk(workdir):
                if File in files:
                    summary_path = root + '/' + File
                    column_name = root.split('/')
                    search_file = open(summary_path, 'r', encoding='utf-8')
                    lines = search_file.readlines()
                    line_clean = []
                    for line in lines:
                        line_cleaned = line.strip('\n')
                        line_clean.append(line_cleaned)
                    line_clean.insert(0, column_name[-1])
                    data.append(line_clean)

            df = pd.DataFrame(data).transpose()
            df.rename(columns=df.iloc[0], inplace=True)
            df = df.iloc[1:]
            df.to_csv('cumulated_variance.csv', sep=' ', quoting=csv.QUOTE_ALL, index=None)
            df.to_csv('cumulated_variance.dat', sep=' ', index=None, header=None)

        return df