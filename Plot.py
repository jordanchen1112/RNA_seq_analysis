import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os 
import cmocean
import seaborn as sns

class RNASeqPlotter:
    def __init__(self, path = './Result_file/', norm_file = 'Normalizatoin_counts.xlsx', deseq2_file = 'DESeq2_result.xlsx', save_path = './Result_img/'):
        self.path = path
        self.norm_file = norm_file
        self.deseq2_file = deseq2_file
        self.save_path = save_path
        self.normc = pd.read_excel(os.path.join(self.path, self.norm_file))
        self.deseq2 = pd.read_excel(os.path.join(self.path, self.deseq2_file))
    
    def is_list(obj):
        list_like_types = (list, np.ndarray, pd.Series, tuple, set)
        return isinstance(obj, list_like_types)
    
    def plot_pca(self, save = False):
        # Set the gene id as index (gene id dont need to use for PCA)
        data_col = list(self.normc.columns)[1:-3]
        data = self.normc[data_col]
        # data.set_index("Gene", inplace = True)
        print(data)

        # Number of principal components to retain      
        n_components = 2  

        # Perform PCA
        pca = PCA(n_components=n_components)
        pca.fit(data.T)
        pca_ = pca.transform(data.T)

        # Labels and color
        labels = data.columns
        color = ['firebrick'] * (len(labels) + 1)
        color[(len(labels) + 1)//2:] = ['cornflowerblue'] * ((len(labels) + 1)//2)
        print(labels)

        # Extract the transformed coordinates for each sample
        x = pca_[:, 0]  # X-axis corresponds to the first principal component
        y = pca_[:, 1]  # Y-axis corresponds to the second principal component

        # Plot the transformed samples
        dpi = 300
        fig = plt.figure(figsize=(8, 8), dpi=dpi)

        plt.scatter(x, y, c = color, s = 150) # c:color/ s:size
        for i, label in enumerate(labels):
            plt.annotate(label, (x[i], y[i]), textcoords="offset points", xytext=(5,10), ha='center',fontsize = 6)

        # Set the plot 
        x_min, x_max = min(x) - 0.5 * abs(min(x)), max(x) + 0.5 * abs(max(x))
        y_min, y_max = min(y) - 0.5 * abs(min(y)), max(y) + 0.5 * abs(max(y))

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCA Analysis')
        
        if save == True:
            plt.savefig(os.path.join(self.save_path, 'PCA.png'), bbox_inches='tight')

        plt.show()

    def plot_volcano(self, p = 0.05, fc = 1.5, interest = None, save = False):
        df = self.deseq2

        # Parameters
        log2fc = df['log2FoldChange']
        neg_log_p_values = -np.log10(df['pvalue'])

        # Set thresholds for significance
        p_value_threshold = -np.log10(p)
        fold_change_threshold = np.log2(fc)

        # Color determination
        colors_map = {'Interest':'orange', 'Upregulated': 'firebrick', 'Downregulated': 'cornflowerblue', 'NS': 'lightgray','Up_N':'rosybrown','Down_N':'lightsteelblue'}
        
        if interest:
            if type(interest) == str:
                Interest = [True if str(interest) in str(desc) else False for desc in df['Description']]
            elif self.is_list(interest):
                Interest = [True if str(desc) in interest else False for desc in df['Orf']]
           
        direction = np.where(log2fc > 0, 'Upregulated', 'Downregulated')
        direction = np.where(Interest, 'Interest', direction)
        direction = np.where((0 < log2fc) & (log2fc < fold_change_threshold), 'Up_N', direction)
        direction = np.where((-fold_change_threshold < log2fc) & (log2fc < 0), 'Down_N', direction)
        direction = np.where(neg_log_p_values > p_value_threshold, direction, 'NS')

        c = [colors_map[d] for d in direction]

        # Plot
        dpi = 300
        fig = plt.figure(figsize=(8, 6), dpi=dpi)

        plt.scatter(log2fc, neg_log_p_values, c = c, alpha=0.9, s = 5)
        plt.axhline(y = p_value_threshold, color='black', linestyle= (0, (5,2.5)),linewidth= 1)
        plt.axvline(x = fold_change_threshold, color='black', linestyle= (0, (5,2.5)),linewidth= 1)
        plt.axvline(x = -fold_change_threshold, color='black', linestyle=(0, (5,2.5)),linewidth = 1 )

        # Add gene name labels for each point

        for i, condition in enumerate(direction):
            if condition == 'Interest':
                plt.text(log2fc[i], neg_log_p_values[i] + 1, df['Gene'][i], fontsize=8, ha='center', va='bottom')

        # Figure settings
 
        plt.xlim(-7, 7)
        plt.ylim(-1, 120)
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-log10(p-value)')
        plt.title(f'Volcano Plot (query = {interest})')

        if save == True:
            plt.savefig(os.path.join(self.save_path, f'Volcano_p_{p}_fc_{fc}.png'), bbox_inches='tight')
        plt.show()

    def plot_ma(self, P_threshold = 0.05, save = False):
        log2avg = np.log2(self.deseq2['baseMean'])
        log2fc = self.deseq2['log2FoldChange']
        p = self.deseq2['pvalue']
        
        colormap = {'Sig':'firebrick','NSig':'dimgrey'}
        color = [colormap[d] for d in np.where(p > P_threshold, 'NSig', 'Sig')]

        dpi = 300
        plt.figure(figsize=(8, 6), dpi=dpi)

        plt.scatter(log2avg, log2fc, s= 1.5, c = color,alpha=1)
        plt.xlim(-5,25)
        plt.ylim(-10,10)
        plt.xlabel('Log2 Average Expression')
        plt.ylabel('Log2 Fold Change')
        plt.title('MA-plot')
        if save == True:
            plt.savefig(os.path.join(self.save_path, f'MA_p_{P_threshold}.png'), bbox_inches='tight')

        plt.show()



    def get_genes(self, query = None, p = 0.05, fc = 2):
        # Get interest genes (query in description / p-vlaue / fold change)
        filter1 = (self.deseq2['log2FoldChange'] > np.log2(fc)) | (self.deseq2['log2FoldChange'] < np.log2(1/fc))
        filter2 = self.deseq2['pvalue'] < p
        filtered_df = self.deseq2[filter1 & filter2]

        name = []
        id = []
        for i in range(len(filtered_df['Description'].values)):
            if query in str(filtered_df['Description'].values[i]):
                name.append(filtered_df['Gene'].values[i])
                id.append(filtered_df['Description'].values[i])

        return name, id
    
    def heat_map(self, query, p = 0.05, fc = 1.5, save = False):
        name, id = self.get_genes(query = query, p = p, fc = fc)
        data_col = list(self.normc.columns)[1:-3]
        data = self.normc[self.normc['Gene'].isin(name)][data_col + ['Gene']]
        data = data[['Gene'] + data_col]
        
        # Drop the 'Gene name' column for the heatmap, but keep it as y-labels
        heatmap_data = data.set_index('Gene')
        heatmap_data = heatmap_data.apply(lambda x: np.log2(x + 1))

        dpi = 300
        num_rows, num_cols = heatmap_data.shape
        plt.figure(figsize=(num_cols, num_rows), dpi = dpi)

        # Store the axis object
        ax = sns.heatmap(heatmap_data, cmap=cmocean.cm.delta, center=0, square=True, yticklabels=True)
        cbar = ax.collections[0].colorbar
        cbar.ax.set_title('Scale Title')
        ax.figure.subplots_adjust(right=0.85)

        plt.title("Gene Expression Heatmap")
        if save == True:
            plt.savefig(os.path.join(self.save_path, f'Heatmap_query_{query}_p_{p}_fc_{fc}.png'), bbox_inches='tight')
    
        plt.show()

if __name__ == '__main__':
    plotter = RNASeqPlotter()
    plotter.plot_pca(save = True)
    plotter.plot_volcano(interest = 'oxidative',save = True )
    plotter.plot_ma(save = True)
    plotter.heat_map(query ='oxidative',save = True)