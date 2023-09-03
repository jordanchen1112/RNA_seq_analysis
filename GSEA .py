import pandas as pd
import numpy as np
import os 
import gseapy as gp
from gseapy.plot import  gseaplot
import matplotlib.pyplot as plt
from gseapy.plot import gseaplot
import tempfile

class GSEA:
    def __init__(self,  file_name = 'DESeq2_result.xlsx', gene_sets = "GO_term_GSEA.gmt", save_path = './Result_file/'):
        self.file_name = os.path.join('./Result_file/' + file_name)
        self.save_path = save_path
        self.gene_sets = os.path.join('./GSEA_CGD_Data/' + gene_sets)
        self.data = pd.read_excel(self.file_name)
        self.gsea_results = None

    def create_temp_gmt_file(self):
        # To solve the ecnding error when reading the gmt file by gseapy
        # Read the file with utf-8 encoding
        with open(self.gene_sets, 'r', encoding='utf-8') as file:
            content = file.read()
            
        # Create a temporary file to save the content
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.gmt')
        with open(temp_file.name, 'w', encoding='utf-8') as file:
            file.write(content)
            
        return temp_file.name 
    
    def Preprocess(self):
        # Remove NA (in padj, Orf)
        data = self.data
        data['padj'] = data['padj'].replace('',pd.NA)
        data['pvalue'] = data['pvalue'].replace('',pd.NA)
        data['Orf'] = data['Orf'].replace('', pd.NA)
        data = data.dropna(subset=['padj', 'pvalue', 'Orf'])

        # Modified orf format (Remove 'orf19.')
        data.loc[:, 'Orf'] = [x[6:] if 'orf' in x else x for x in data['Orf']]
        
        return data

    def Prerank_analysis(self):
        fillter_data = self.Preprocess()
   
        # Creat Rank file for GSEA (Pre rank)
        fillter_data['Rank'] = fillter_data['log2FoldChange'] # Create 'Rank' column in the dataframe 
        fillter_data = fillter_data.sort_values('Rank', ascending = False) # 排序   
        ranking = fillter_data[['Orf','Rank']]
        print(ranking)

        # Gene sets
        temp_gene_sets = self.create_temp_gmt_file()

        # 執行 GSEA Prerank (Enrichment)
        pre_results = gp.prerank(
            rnk=ranking,
            gene_sets= temp_gene_sets,
            seed = 42,
            verbose=True,
        )
        self.gsea_results = pre_results

    def Prerank_result(self):
        result = self.gsea_results
        out = []
        for t in result.results:
            p = result.results[t]['pval']
            fdr = result.results[t]['fdr']
            nes = result.results[t]['nes']
            es = result.results[t]['es']
            lead_genes = result.results[t]['lead_genes']
            matched_genes = result.results[t]['matched_genes']
            gene_p = result.results[t]['gene %']
            out.append([t, p,fdr,nes,es,abs(nes),lead_genes,matched_genes,gene_p])

            prerank_df = pd.DataFrame(out,columns =['Term','p_value','fdr','nes','es',
                                                    'abs','lead_genes','matched_genes','gene%'])
       
        prerank_df.to_excel(os.path.join(self.save_path, f'GSEA_Prerank.xlsx'), index=False)

        return prerank_df

    def EnrichmentPlot(self, term = None, save_path = None):
        if not save_path:
            save_path = self.save_path
        if term == None:
            print('Please Enter the GO term you want to plot.')

        result = self.gsea_results
        all_gene_sets = result.results.keys()
        terms = [s for s in all_gene_sets if term in s]
        for term in terms:
            print(f'ploting {term}')
            if save_path != None:
                gseaplot(rank_metric=result.ranking, term=terms[0], ofname= os.path.join(save_path, f'{term}.pdf'), **result.results[terms[0]])
            else:
                gseaplot(rank_metric=result.ranking, term=terms[0], ofname=f'{term}.pdf', **result.results[terms[0]])

if __name__ == '__main__':
    Gsea = GSEA()
    Gsea.Preprocess()
    Gsea.Prerank_analysis()
    Gsea.Prerank_result()
    Gsea.EnrichmentPlot(term = 'oxid', save_path = './Result_img/')

