import numpy as np
import pandas as pd
import os 
import plotly.express as px
import plotly.graph_objects as go

class InteractivePlotter:
    def __init__(self, path='./Result_file/', norm_file='Normalizatoin_counts.xlsx', deseq2_file='DESeq2_result.xlsx', save_path='./Result_img/'):
        self.path = path
        self.norm_file = norm_file
        self.deseq2_file = deseq2_file
        self.save_path = save_path
        self.normc = pd.read_excel(os.path.join(self.path, self.norm_file))
        self.deseq2 = pd.read_excel(os.path.join(self.path, self.deseq2_file))

    def wrap_text(self, text, max_length=50):
        """
        Wraps the input text every max_length characters.
        """
        text = str(text)
        words = text.split()
        lines = []
        current_line = []

        for word in words:
            # If adding the new word doesn't exceed the maximum line length, add the word to the current line
            if len(' '.join(current_line) + ' ' + word) <= max_length:
                current_line.append(word)
            else:
                # Otherwise, start a new line
                lines.append(' '.join(current_line))
                current_line = [word]
        
        # Add the last line
        lines.append(' '.join(current_line))
        
        return '<br>'.join(lines)

    def is_list(self, obj):
        list_like_types = (list, np.ndarray, pd.Series, tuple, set)
        return isinstance(obj, list_like_types) 

    def interactive_volcano(self, p=0.05, fc=1.5, interest=None):
        df = self.deseq2

        log2fc = df['log2FoldChange']
        neg_log_p_values = -np.log10(df['pvalue'])
        p_value_threshold = -np.log10(p)
        fold_change_threshold = np.log2(fc)
        df['WrappedDescription'] = df['Description'].apply(self.wrap_text)
        print(df['WrappedDescription'][11])
        colors_map = {'Interest':'orange', 'Upregulated': 'firebrick', 'Downregulated': 'cornflowerblue', 'NS': 'lightgray','Up_N':'rosybrown','Down_N':'lightsteelblue'}
        
        Interest = [False]*len(df)
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

        df['color'] = [colors_map[d] for d in direction]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
                x=log2fc,
                y=neg_log_p_values,
                mode='markers',
                marker=dict(
                    color=df['color'],
                    size=6
                ),
                text = 'Gene name: ' + df['Gene'] + '<br>' + 'log2FC: ' + round(df['log2FoldChange'], 2).astype(str) 
                + '<br>' + 'P-value: ' + round(df['padj'], 6).astype(str)
                + '<br>' + df['WrappedDescription'],
                hoverinfo='text'
            ))

        # Add lines for thresholds
        fig.add_shape(type="line", x0=log2fc.min(), x1=log2fc.max(), y0=p_value_threshold, y1=p_value_threshold, line=dict(color="Black", width=2, dash="dash"))
        fig.add_shape(type="line", x0=fold_change_threshold, x1=fold_change_threshold, y0=0, y1=max(neg_log_p_values), line=dict(color="Black", width=2, dash="dash"))
        fig.add_shape(type="line", x0=-fold_change_threshold, x1=-fold_change_threshold, y0=0, y1=max(neg_log_p_values), line=dict(color="Black", width=2, dash="dash"))

        # Adjust layout if needed
        fig.update_layout(
            width=800,
            height = 600,  
            xaxis_title='Log2 Fold Change',
            yaxis_title='-log10(p-value)',
            xaxis_range=[-8, 8],  
            yaxis_range=[-15, 120],
            uirevision='constant'
        )
        
        fig.show()
 
        # Display the plot
        return fig

if __name__ == '__main__':
    test = InteractivePlotter()
    test.interactive_volcano(interest = 'oxidative')


