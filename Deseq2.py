# Differential_expression.py

import pandas as pd
import numpy as np
import os 
import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2 import preprocessing
import re
from Translator import Translator
from Preprocessing import Preprocessing

class DifferentialExpression:
    def __init__(self, output_path, processed_data):
        # Read preprocessed count file
        self.output_path = output_path
        self.processed_data = processed_data
        self.counts = None
        self.matadata= None
        self.estimated_gene = None  # Estimated gene expression dataframe
        self.result = None  # DESeq2 result dataframe

    def prepare_count_data(self):
        # Convert your preprocessed data into the format required by pydeseq2
        # Get Processed Read counts
        counts = self.processed_data[['transcript_id', 'Rap1Del_2', 'Rap1Del_3', 'Rap1Del_4', 'Rap1WT_1','Rap1WT_2', 'Rap1WT_3']] 
        counts.set_index("transcript_id", inplace = True)
        counts = counts.astype(int) # Convert values to "int"
        counts = counts.T # 將 Dataframe 轉置
        self.counts = counts
        
    def prepare_meta_data(self):
        # <metadata> contains the condition design in this experement
        metadata = pd.DataFrame(zip(self.counts.index, ['Rap1Del', 'Rap1Del', 'Rap1Del', 'Rap1WT', 'Rap1WT', 'Rap1WT']), columns = ['sample','condition'] )
        metadata.set_index('sample', inplace =True)
        self.metadata = metadata

    def run_deseq2_analysis(self):
        dds = DeseqDataSet(
            counts=self.counts,
            metadata =self.metadata,
            design_factors="condition",
            refit_cooks=True,
            )
        
        # Run DDS
        dds.deseq2()
        deseq2_result = DeseqStats(dds, contrast = ('condition','Rap1Del','Rap1WT'), alpha=0.05, cooks_filter=True, independent_filter=True)
        deseq2_result.summary()

        # Get the DESeq2 result
        result_df = deseq2_result.results_df
        
        # Get the ORF, Gene name, Discription by pre-build dictionary  
        id_orf_transform, id_gene_name_transform, id_discription_transform = Translator()
        result_df.insert(0, 'ID', result_df.index)
        result_df.insert(1, 'Orf', result_df['ID'].map(id_orf_transform))
        result_df.insert(2, 'Gene', result_df['ID'].map(id_gene_name_transform))
        result_df.insert(3, 'Discription', result_df['ID'].map(id_discription_transform))
        self.result = result_df

        # Get the Estimated Gene Expression
        estimated_gene_expression = dds.layers['normed_counts']
        estimated_df = pd.DataFrame(estimated_gene_expression.T, 
                                                    index=dds.var.index, columns=dds.obs.index) # Turn to dataframe
        estimated_df.insert(0, 'ID', estimated_df.index)
        estimated_df.insert(1, 'Orf', estimated_df['ID'].map(id_orf_transform))
        estimated_df.insert(2, 'Gene', estimated_df['ID'].map(id_gene_name_transform))
        self.estimated_gene = estimated_df

    def save_results(self):
        self.result.to_excel(os.path.join(self.output_path, 'DESeq2_result.xlsx'), index=False)
        self.estimated_gene.to_excel(os.path.join(self.output_path, 'Estimated_gene.xlsx'), index=False)

if __name__ == '__main__':
    path = './Rap1/Data/'
    file_name = 'Preprocess_TS_Count.xlsx'
    processed_data = pd.read_excel(os.path.join(path,file_name))
    DE = DifferentialExpression(output_path = path, processed_data = processed_data)
    DE.prepare_count_data()
    DE.prepare_meta_data()
    DE.run_deseq2_analysis()
    DE.save_results()
