import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import os
#from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

class Plot:

    def bar_plot(self, xlabel, ylabel, zlabel, xxlabels, yylabels, view, save, variance=False, cumulated_variance=False):

        if variance == True:
            df = pd.read_csv('variance.dat', sep=' ', names=list(range(15)))

            self.threed_bar_plot(df, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, xxlabels=xxlabels, yylabels=yylabels, view=True, save=False)

        if cumulated_variance == True:
            df = pd.read_csv('cumulated_variance.dat', sep=' ')
            print(os.getcwd())
            print(len(df.columns))

            self.threed_bar_plot(df, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, xxlabels=xxlabels, yylabels=yylabels, view=True, save=False)
            
        return

    def threed_bar_plot(self, filename, xlabel, ylabel, zlabel, xxlabels, yylabels, view=True, save=False): #######
        
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(111, projection='3d')

        data = filename.values.tolist()
        x = []
        y = []
        z = []

        for i in range(len(data)):
            for j in range(len(data[i])):
                x.append(i)
                y.append(j)
                z.append(data[i][j])

        top = x + y
        bottom = np.zeros_like(x)
        width = depth = 0.6

        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        ax1.set_zlabel(zlabel)
        xlabels = np.array(xxlabels)
        xpos = np.arange(xlabels.shape[0])
        ylabels = np.array(yylabels)
        ypos = np.arange(ylabels.shape[0])
        zpos = data

        dx=0.5
        dy=0.5
        dz=zpos

        ax1.w_xaxis.set_ticks(xpos + dx/2.)
        ax1.w_xaxis.set_ticklabels(xlabels)

        ax1.w_yaxis.set_ticks(ypos + dy/2.)
        ax1.w_yaxis.set_ticklabels(ylabels)

        cmap = cm.get_cmap('hsv')
        max_height = np.max(z)  
        min_height = np.min(z)
        rgba = [cmap((k-min_height)/max_height) for k in z] 
        norm = cm.colors.Normalize(vmin=min_height, vmax=max_height, clip=False)

        ax1.bar3d(x, y, bottom, width, depth, z, shade=True, edgecolor='black', linewidth=0.3, color=rgba)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))

        ax1.set_title(' ')

        if view == True: 
            print('Now you are visualizing the plot...')
            plt.show()
        if save == True:
            print('Your plot has been saved.')
            plt.savefig('plot_3d.png')

        return

    def explained_variance_plot(self, explained_variance):

        labels, ys = zip(*explained_variance)
        xs = np.arange(len(labels)) 
        width = 0.4

        plt.bar(xs, ys, width, align='center')

        plt.xticks(xs, labels) #Replace default x-ticks with xs, then replace xs with labels
        plt.yticks(ys)

        plt.xlabel("pc")
        plt.ylabel("Explained Variance")
        plt.title("Explained Variance")

        plt.savefig('explained_variance.png')

        return

    def variance_plot(self, variance, pc, zeroth_dimension, nth_dimension):

        plt.plot(pc.variance[zeroth_dimension:nth_dimension])
        plt.xlabel('Principal component')
        plt.ylabel('Variance');

        plt.savefig('variance.png')

        return

    def cumulated_variance_plot(self, cumulated_variance, pc, zeroth_dimension, nth_dimension):

        plt.plot(pc.cumulated_variance[zeroth_dimension:nth_dimension])
        plt.xlabel('Principal component')
        plt.ylabel('Cumulated variance');

        plt.savefig('cumulated_variance.png')

        return

    def time_plot(self, df):
        
        g = sns.PairGrid(df, hue='Time (ps)',
                         palette=sns.color_palette('Oranges_d',
                                                   n_colors=len(df)))
        g.map(plt.scatter, marker='.')

        g.savefig('pairGrid.png')

        return

########## UNDER DEVELOPMENT #############
    # def time_scatterplot(self, df):

    #     g = sns.scatterplot(df)

    #     g.savefig('pcGraph.png') 
        
    #     return