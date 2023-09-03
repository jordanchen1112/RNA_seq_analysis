# Use the FASTA file downloaded from CGD to link the "transcript_id" with [ORF/gene_name/decription]

import re
import pandas as pd
import os

def Translator(fasta_file='./GSEA_CGD_Data/C_albicans_SC5314_A22_current_orf_coding.fasta'):

    id_orf_transform = {}
    id_gene_name_transform = {}
    id_description_transform = {}
    gene_info = {}

    # Read FASTA file
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespaces
            if line.startswith(">"):  # Description starts with ">" in FASTA file
                # Use RE to catch the ['ID','Gene','Orf']
                match = re.match(r'>(\w+_\w+_\w+) (\w+).*\((orf[\d\.]+)\)', line)
                if match:
                    transcript_id = match.group(1)
                    gene_name = match.group(2)
                    orf = match.group(3)
                
                if transcript_id[-1] == 'B':
                    continue
                else:
                    transcript_id = transcript_id[:-2]
                    
                if re.match(r'C[0-9R]_\d+[CW]_[AB]', gene_name):
                    gene_name = orf
                
                if not gene_name:
                    gene_name = orf

                id_orf_transform[transcript_id] = orf
                id_orf_transform[orf] = transcript_id
                id_gene_name_transform[transcript_id] = gene_name
                id_gene_name_transform[gene_name] = transcript_id

                gene_description = line.split(";", 1)[1].strip()
                if gene_description:
                    if '(orf' in gene_description:
                        gene_description = gene_description.split(")", 1)[1].strip()
                else:
                    gene_description = ''

                id_description_transform[transcript_id] = gene_description
                id_description_transform[gene_description] = transcript_id

                gene_info[transcript_id] = [transcript_id, gene_name, orf, gene_description]

    gene_info_df = pd.DataFrame.from_dict(gene_info, orient='index', columns=['ID', 'Gene', 'Orf', 'Description'])
    gene_info_df.to_excel(os.path.join('./GSEA_CGD_Data/', 'gene_info.xlsx'), index=False)

    return gene_info_df, id_orf_transform, id_gene_name_transform, id_description_transform

if __name__ == "__main__":
    
    gene_info_df, orf_transform, gene_name_transform, description_transform = Translator()
    print(gene_info_df)
    
    # print(orf_transform)
    # print(gene_name_transform)
    # print(description_transform)
    
    # # Get the ORF, Gene name, Discription by pre-build dictionary  
    # id_orf_transform, id_gene_name_transform, id_discription_transform = Translator()
    # result_df.insert(0, 'ID', result_df.index)
    # result_df.insert(1, 'Orf', result_df['ID'].map(id_orf_transform))
    # result_df.insert(2, 'Gene', result_df['ID'].map(id_gene_name_transform))
    # result_df.insert(3, 'Discription', result_df['ID'].map(id_discription_transform))
    # self.result = result_df

    # Get the Estimated Gene Expression
    # estimated_gene_expression = dds.layers['normed_counts']
    # estimated_df = pd.DataFrame(estimated_gene_expression.T, 
    #                                             index=dds.var.index, columns=dds.obs.index) # Turn to dataframe
    # estimated_df.insert(0, 'ID', estimated_df.index)
    # estimated_df.insert(1, 'Orf', estimated_df['ID'].map(id_orf_transform))
    # estimated_df.insert(2, 'Gene', estimated_df['ID'].map(id_gene_name_transform))
    # self.estimated_gene = estimated_df