import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os 
from Preprocessing import Preprocessing

class RNASeqPlotter:
    def __init__(self, path = './Rap1/Data/'):
        self.path = path
        self.data = pd.read_excel(os.path.join(self.path, 'DESeq2_result.xlsx'))
        self.estimated_gene = pd.read_excel(os.path.join(self.path, 'Preprocess_TS_TPM.xlsx'))

    def is_list(obj):
        list_like_types = (list, np.ndarray, pd.Series, tuple, set)
        return isinstance(obj, list_like_types)
    
    def plot_pca(self, save = False):
        norm_counts = self.estimated_gene[['Gene ID', 'Rap1Del_2', 'Rap1Del_3', 'Rap1Del_4', 'Rap1WT_1', 'Rap1WT_2', 'Rap1WT_3']]
        norm_counts.set_index("Gene ID", inplace = True)  

        # Number of principal components to retain      
        n_components = 2  

        # Perform PCA
        pca = PCA(n_components=n_components)
        pca.fit(norm_counts.T)
        pca_ = pca.transform(norm_counts.T)

        # Labels and color
        labels = ['Rap1Del_2', 'Rap1Del_3', 'Rap1Del_4', 'Rap1WT_1','Rap1WT_2', 'Rap1WT_3']
        color = ['firebrick','firebrick','firebrick',
                 'cornflowerblue','cornflowerblue','cornflowerblue']

        # Extract the transformed coordinates for each sample
        x = pca_[:, 0]  # X-axis corresponds to the first principal component
        y = pca_[:, 1]  # Y-axis corresponds to the second principal component

        # Plot the transformed samples
        dpi = 300
        fig = plt.figure(figsize=(8, 8), dpi=dpi)

        plt.scatter(x, y, c = color, s = 30)
        for i, label in enumerate(labels):
            plt.annotate(label, (x[i], y[i]), textcoords="offset points", xytext=(5,10), ha='center',fontsize = 6)

        plt.ylim(-40000, 40000)
        plt.xlim(-40000, 40000)

        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA Analysis')
        plt.show()

        if save == True:
            plt.savefig(os.path.join(self.path, 'PCA.png'), bbox_inches='tight')

    def plot_volcano(self, p = 0.05, fc = 1.5, label = 'interest', interest = None, save = False):
 
        df = self.data    
        # Parameters
        log2fc = df['log2FoldChange']
        neg_log_p_values = -np.log10(df['pvalue'])

        names = np.where(df['Gene'].notnull(), df['Gene'], df['Orf'])

        # Set thresholds for significance
        p_value_threshold = -np.log10(p)
        fold_change_threshold = np.log2(fc)

        # Color determination
        direction = np.where(log2fc > 0, 'Upregulated', 'Downregulated')
        direction = np.where((0 < log2fc) & (log2fc < fold_change_threshold), 'Up_N', direction)
        direction = np.where((-fold_change_threshold < log2fc) & (log2fc < 0), 'Down_N', direction)
        direction = np.where(neg_log_p_values > p_value_threshold, direction, 'NS')

        colors_map = {'Upregulated': 'firebrick', 'Downregulated': 'cornflowerblue', 'NS': 'lightgray','Up_N':'rosybrown','Down_N':'lightsteelblue'}
        if interest:
            if type(interest) == str:
                Interest = [True if str(interest) in str(desc) else False for desc in df['Discription']]
            elif self.is_list(interest):
                Interest = [True if str(desc) in interest else False for desc in df['Orf']]
                direction = np.where(Interest, 'Interest', direction)
                colors_map['Interest'] = 'orange'
        c = [colors_map[d] for d in direction]

        # Plot
        dpi = 300
        fig = plt.figure(figsize=(8, 6), dpi=dpi)

        plt.scatter(log2fc, neg_log_p_values, c = c, alpha=0.9, s = 5)
        plt.axhline(y = p_value_threshold, color='black', linestyle= (0, (5,2.5)),linewidth= 1)
        plt.axvline(x = fold_change_threshold, color='black', linestyle= (0, (5,2.5)),linewidth= 1)
        plt.axvline(x = -fold_change_threshold, color='black', linestyle=(0, (5,2.5)),linewidth = 1 )

        # Add gene name labels for each point
        if label == 'all':
            for i, gene_name in enumerate(names):
                if 0 < log2fc[i] < 2 and neg_log_p_values[i] < 25:
                    continue
                if 0 > log2fc[i] > - 2 and neg_log_p_values[i] < 25:
                    continue
                if neg_log_p_values[i] < 5:
                    continue
                plt.text(log2fc[i], neg_log_p_values[i] + 1, gene_name, fontsize=5, ha='center', va='bottom')

        elif label == 'interest':
            if interest:
                for i, condition in enumerate(direction):
                    if neg_log_p_values[i] < p_value_threshold:
                        continue
                    elif 0 < log2fc[i] < fold_change_threshold or -fold_change_threshold < log2fc[i] < 0:
                        continue
                    elif condition == 'Interest':
                        plt.text(log2fc[i], neg_log_p_values[i] + 1, names[i], fontsize=5, ha='center', va='bottom')

        # Figure settings
        plt.xlim(-7, 7)
        plt.ylim(-1, 120)
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        plt.title('Volcano Plot')
        plt.show()

        if save == True:
            plt.savefig(os.path.join(self.path, 'Volcano.png'), bbox_inches='tight')

    def plot_ma(self, P_threshold = 0.05):
        log2avg = np.log2(self.data['baseMean'])
        log2fc = self.data['log2FoldChange']
        p = self.data['pvalue']
        
        colormap = {'Sig':'firebrick','NSig':'dimgrey'}
        color = [colormap[d] for d in np.where(p > P_threshold, 'NSig', 'Sig')]

        dpi = 300
        fig = plt.figure(figsize=(8, 6), dpi=dpi)

        plt.scatter(log2avg, log2fc, s= 1.5, c = color,alpha=1)
        plt.xlim(-5,25)
        plt.ylim(-10,10)
        plt.xlabel('Log2 Average Expression')
        plt.ylabel('Log2 Fold Change')
        plt.title('MA-plot')
        plt.show()

if __name__ == '__main__':
    plotter = RNASeqPlotter()
    plotter.plot_pca(save = True)