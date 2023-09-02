import pandas as pd
import numpy as np
import os 
import gseapy as gp
from gseapy.plot import  gseaplot
import matplotlib.pyplot as plt
from gseapy.plot import gseaplot

class GSEA:
    def __init__(self, path, file_name, gene_sets = "./Rap1/GSEA/GO_term_GSEA.gmt"):
        # Read DESeq2_result file
        self.path = path
        self.file_name = file_name + '.xlsx'
        self.gene_sets = gene_sets
        self.data = pd.read_excel(os.path.join(self.path, self.file_name))
        self.gsea_results = None
    
    def Preprocess(self):
        # Remove NA (in padj, Orf)
        data = self.data
        data['padj'] = data['padj'].replace('',pd.NA)
        data['Orf'] = data['Orf'].replace('', pd.NA)
        data = data.dropna(subset=['padj', 'Orf'])

        # Modified orf format (Remove 'orf19.')
        data['Orf'] = [x[6:] if 'orf' in x else x for x in data['Orf']]
        
        return data

    def Prerank_analysis(self):
        data = self.Preprocess()
   
        # Creat Rank file for GSEA (Pre rank)
        data['Rank'] = data['log2FoldChange'] # Create 'Rank' column in the dataframe 
        data = data.sort_values('Rank', ascending = False) # 排序   
        ranking = data[['Orf','Rank']]
    
        # Gene set file
        gene_sets = self.gene_sets

        # 執行 GSEA Prerank (Enrichment)
        results = gp.prerank(
            rnk=ranking,
            gene_sets= gene_sets,
            seed = 42,
            verbose=True,
        )
        self.gsea_results = results

    def Prerank_result(self , save_path = None):
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
        if save_path != None:
            prerank_df.to_excel(os.path.join(save_path, 'GSEA_Prerank_result.xlsx'), index=False)

        return prerank_df

    def EnrichmentPlot(self, term = None, save_path = None):
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
    Gsea = GSEA(path = './Rap1/Data/', file_name = 'DESeq2_result', gene_sets = "./Rap1/GSEA/GO_term_GSEA.gmt")
    Gsea.Preprocess()
    Gsea.Prerank_analysis()
    Gsea.Prerank_result(save_path = './Rap1/GSEA/')
    Gsea.EnrichmentPlot(term = 'oxid', save_path = './Rap1/GSEA/')

