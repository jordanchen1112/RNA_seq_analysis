# RNA_seq_analysis
- This project provides tools and scripts for the analysis of RNA Sequencing data. It covers various aspects of the analysis pipeline, from preprocessing raw data to performing differential expression analysis and gene set enrichment analysis (GSEA).

# Files Description
1. Deseq2.py: Contains the DifferentialExpression class which performs differential expression analysis using the DESeq2 approach.

2. GSEA .py: Contains the GSEA class which performs Gene Set Enrichment Analysis and creates relevant plots.

3. Plot.py: Contains the RNASeqPlotter class which offers various plotting functions to visualize RNA Sequencing analysis results.

4. Preprocessing.py: Contains the Preprocessing class which handles the preprocessing of raw RNA sequencing data.

5. Translator.py: Provides the Translator function which links "transcript_id" with attributes like ORF, gene name, and description using a FASTA file.

## RNA Seq Analysis: Step-by-Step Guide

### 1. Setup and Installation
- Clone the repository to your local machine.
    ```sh
    git clone <repository_url>
    ```
- If git was not installed, go to <> code/Donwload Zip.<br> 
    This could download all the files.

- Install the necessary Python libraries:
    ```sh
    pip install -r requirements.txt
    ```

### 2. Preprocessing Raw Data (Preprocessed.py)
- **Input**: Raw RNA sequencing data (`your_raw_data_file.xlsx`)<br>
    - Put the raw data files (Read counts) in the folder 'Raw_data' and the columns should the same as that in 'Example_file.xlsx'<br>

- **Output**: Preprocessed data (`Preprocess_ + {your_raw_data_file}.xlxs`)
    - Output excel file will default save in the ./Result_file/ 
    - This output excel file could be further used for Deseq2.py

    ```python
    from Preprocessing import Preprocessing
    
    # Initialize and preprocess
    preprocess = Preprocessing(file_name = 'your_raw_data_file.xlsx')
    preprocess.read_file()
    preprocess.average_alleles()  # The proccessing processes are optional
    preprocess.save_file()
    ```

### 3. Differential Expression Analysis (Deseq2.py)
- **Input**: Preprocessed data (`Preprocess_result.xlsx`) 
    - Get from **Preprocessed.py**
- **Output**: 
    1. Differential expression results (`DESeq2_result.xlsx`)
    2. Normalized counts results (`Normalizatoin_counts.xlsx`)
    - Output excel file will default save in the ./Result_file/ <br>

    ```python
    from Deseq2 import DifferentialExpression
    
    # Initialize and perform analysis
    DE = DifferentialExpression(file_name = 'Preprocess_TS_Count.xlsx')
    DE.prepare_dds()
    DE.run_deseq2_analysis()
    DE.add_info()
    DE.save_results()

    ```

### 4. Gene Set Enrichment Analysis (GSEA.py)
- **Input**: 
    1. Differential expression results (`DESeq2_result.xlsx`)
    2. Gene sets file (.gmt)
    - `DESeq2_result.xlsx` could get from **Deseq2.py**

- **Output**: GSEA results (`gsea_results_file.ext`)

    ```python
    from GSEA import GSEA
    
    # Initialize and perform GSEA
    Gsea = GSEA(file_name = 'DESeq2_result.xlsx', gene_sets = "GO_term_GSEA.gmt")
    Gsea.Preprocess()
    Gsea.Prerank_analysis()
    Gsea.Prerank_result()
    Gsea.EnrichmentPlot(term = 'oxid', save_path = './Result_img/')


### 5. Data Visualization
- **Input**:     
    1. Differential expression results (`DESeq2_result.xlsx`)
    2. Normalized counts results (`Normalizatoin_counts.xlsx`) for PCA
    - Output excel file will default save in the ./Result_img/ <br>
- **Output**: Various plots (e.g., heatmap, PCA)

    ```python
    from Plot import RNASeqPlotter
    
    # Initialize plotter
    plotter = RNASeqPlotter(norm_file = 'Normalizatoin_counts.xlsx', deseq2_file = 'DESeq2_result.xlsx')
    
    # Create plots
    plotter.plot_pca(save = True)
    plotter.plot_volcano(interest = 'oxidative',save = True )
    plotter.plot_ma(save = True)
    plotter.heat_map(query ='oxidative',save = True)

### 6. Data Translation (Optional)
- **Input**: 
  - FASTA file (`your_fasta_file.fasta`)
- **Output**: Translated data file with linked attributes (`translated_data_file.ext`)

    ```python
    from Translator import Translator
    
    # Translate data
    translated_data = Translator('path_to_your_fasta_file.fasta', 'path_to_your_data_file.ext')
    ```

### 7. Review and Analysis
- Go through the generated results, plots, and any translated data to derive insights and conclusions.
