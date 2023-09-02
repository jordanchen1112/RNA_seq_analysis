# Preprocessing.py

import pandas as pd
import os 

class Preprocessing:
    def __init__(self, path, file_name):
        self.path = path
        self.file_name = file_name
        self.data = pd.read_excel(os.path.join(path, file_name))
    
    def average_alleles(self, save = False):
        # Remove the 'A' & 'B' of the transcription_id
        self.data['transcript_id'] = self.data['transcript_id'].str[:-4]
        cols = ['Rap1Del_2', 'Rap1Del_3', 'Rap1Del_4', 'Rap1WT_1', 'Rap1WT_2', 'Rap1WT_3']
        merged_df = self.data.groupby('transcript_id').agg({**{col: 'mean' for col in cols}}).reset_index()
        self.data = merged_df

        if save == True:
            name = 'Preprocess_' + self.file_name
            self.data.to_excel(os.path.join(self.path, name), index=False)

if __name__ == '__main__':
    Preprocessor = Preprocessing(path = './Rap1/Data/', file_name = 'TS_Count.xlsx')
    Preprocessor.average_alleles(save =True)
