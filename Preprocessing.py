# Preprocessing.py

import pandas as pd
import os 

class Preprocessing:
    def __init__(self, path = './Raw_data/', file_name = 'TS_Count.xlsx', save_path = './Result_file/'):
        self.path = path
        self.file_name = file_name
        self.save_path = save_path
        self.id = None
        self.sample_cols = None
        
    def read_file(self):
        self.data = pd.read_excel(os.path.join(self.path, self.file_name))
        cols = list(self.data.columns)
        self.id = str(cols[0])
        self.sample_cols = list(cols[1:])
       
    def average_alleles(self):
        # Remove the 'A' & 'B' of the transcription_id
        self.data[self.id] = self.data[self.id].str[:-4]
        self.data = self.data.groupby(self.id).agg({**{col: 'mean' for col in self.sample_cols}}).reset_index()
 
    def remove_zero(self):
        self.data = self.data[(self.data != 0).all(axis=1)]

    def save_file(self):
        save_name = 'Preprocess_' + self.file_name
        self.data.to_excel(os.path.join(self.save_path, save_name), index=False)

if __name__ == '__main__':
    Rap1 = Preprocessing(file_name = 'TS_Count.xlsx')
    Rap1.read_file()
    Rap1.average_alleles()
    # Rap1.remove_zero()
    Rap1.save_file()
