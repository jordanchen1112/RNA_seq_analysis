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
    def __init__(self, path = './Result_file/', file_name = 'Preprocess_TS_Count.xlsx', save_path = './Result_file/'):
        # Read preprocessed count file
        self.path = path
        self.file_name = file_name
        self.save_path = save_path
        self.pre_data = pd.read_excel(os.path.join(self.path, self.file_name))
        self.counts = None
        self.matadata= None
        self.norm_counts = None  # Normalization read counts 
        self.result = None  # DESeq2 result 

    def prepare_dds(self):
        # Prepare the files (counts / meta data) required by pydeseq2
        
        # Get Columns Names from Processed Excel File 
        cols = list(self.pre_data.columns)
        self.id = str(cols[0])
        self.sample_cols = list(cols[1:])

        # <Count> data contains the gene id and original read counts
        counts = self.pre_data[cols] 
        counts.set_index(f'{self.id}', inplace = True)
        counts = counts.astype(int) # Convert values to "int"
        counts = counts.T # 將 Dataframe 轉置
        self.counts = counts
        
        # <metadata> contains the condition design in this experement
        sample = self.sample_cols
        condition = [x[:-2] for x in self.sample_cols]
        metadata = pd.DataFrame(zip(sample, condition ), 
                                columns = ['sample','condition'] )
        metadata.set_index('sample', inplace =True)
        self.metadata = metadata
        print(self.metadata)

    def run_deseq2_analysis(self):
        dds = DeseqDataSet(
            counts=self.counts,
            metadata =self.metadata,
            design_factors="condition",
            refit_cooks=True,
            )
        
        # Run DDS
        dds.deseq2()
        condition = [x[:-2] for x in self.sample_cols]
        conditionA = condition[0]
        conditoinB =condition[-1]
        deseq2_result = DeseqStats(dds, contrast = ('condition',str(conditionA),str(conditoinB)),
                                   alpha=0.05, cooks_filter=True, independent_filter=True)
        # Run statistical result
        deseq2_result.summary()

        # Get the DESeq2 result
        self.result = deseq2_result.results_df

        # Get Normalization counts and Turn into dataframe
        self.norm_counts = pd.DataFrame(dds.layers['normed_counts'].T, index=dds.var.index, columns=dds.obs.index) 
    
    def add_info(self):
        fasta_file='./GSEA_CGD_Data/C_albicans_SC5314_A22_current_orf_coding.fasta'
        gene_info_df, orf_transform, gene_name_transform, description_transform = Translator(fasta_file = fasta_file)
        self.result = self.result.merge(gene_info_df, on='ID', how='left')
        self.norm_counts = self.norm_counts.merge(gene_info_df, on='ID', how='left')

    def save_results(self):
        self.result.to_excel(os.path.join(self.save_path, 'DESeq2_result.xlsx'), index = False)
        self.norm_counts.to_excel(os.path.join(self.save_path, 'Normalizatoin_counts.xlsx'), index = False)

if __name__ == '__main__':
    DE = DifferentialExpression()
    DE.prepare_dds()
    DE.run_deseq2_analysis()
    DE.add_info()
    DE.save_results()
