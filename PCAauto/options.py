import pandas as pd
import numpy as np
import MDAnalysis as mda

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
                          columns=['PC{}'.format(i+1) for i in range(n_dim)])
        df['Time (ps)'] = df.index * self.u.trajectory.dt
        df.head()

        df.to_csv('dataFrame.csv', index=True)
        df.to_excel('dataFrame.xlsx', engine='xlsxwriter') 

        return df, transformed

    def visualize_movement(self, transformed): #This is bugged

        pc1 = self.pc.p_components[:, 0]
        trans1 = transformed[:, 0]
        projected = np.outer(trans1, pc1) + self.pc.mean
        coordinates = projected.reshape(len(trans1), -1, 3)

        proj1 = mda.Merge(backbone)
        proj1.load_new(coordinates, order="fac")

        proj1.atoms.write('pca1.pdb')
        proj1.atoms.write('pca1.dcd', frames='all')

        return proj1