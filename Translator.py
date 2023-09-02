# Use the FASTA file downloaded from CGD to link the "transcript_id" with [ORF/gene_name/decription]

import re

def Translator(fasta_file = './Rap1/GSEA/C_albicans_SC5314_A22_current_orf_coding.fasta'):

    id_orf_transform = {}
    id_gene_name_transform = {}
    id_discription_transform = {}

    # Read FASTA file
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip() # Remove leading/trailing whitespaces
            if line.startswith(">"): # Dicription starts with ">" in FASTA file
                match = re.match(r'>(\w+)_\w+ (\w+).*\((orf[\d\.]+)\)', line) # use RE to catch the [orf, gene name, transcript id]
                if match:
                    transcript_id = match.group(1)
                    orf = match.group(3)
                    gene_name = match.group(2)
                    id_orf_transform[transcript_id] = orf
                    id_orf_transform[orf] = transcript_id

                    if re.match(r'C[0-9R]_\d+[CW]_[AB]', gene_name):
                        id_gene_name_transform[transcript_id] = ''
                    else:
                        id_gene_name_transform[transcript_id] = gene_name
                        id_gene_name_transform[gene_name] = transcript_id
                else:
                    id_orf_transform[transcript_id] = ''
                    id_gene_name_transform[transcript_id] = ''

                gene_description = line.split(";",1)[1].strip()

            if gene_description:
                if '(orf' in gene_description:
                    gene_description = gene_description.split(")", 1)[1].strip()
                else:
                    gene_description = gene_description
                id_discription_transform[transcript_id] = gene_description
                id_discription_transform[gene_description] = transcript_id
            else:
                id_discription_transform[transcript_id] = ''

    return id_orf_transform, id_gene_name_transform, id_discription_transform
    
if __name__ == "__main__":
    fasta_file = './Rap1/GSEA/C_albicans_SC5314_A22_current_orf_coding.fasta'
    orf_transform, gene_name_transform, description_transform = Translator(fasta_file)
    # print(orf_transform)
    # print(gene_name_transform)
    # print(description_transform)