# scRNAseq Seurat Analysis Shiny App
This Shiny app facilitates the analysis of single-cell RNA sequencing (scRNAseq) data using the Seurat package in R.

##  Features

 Key features include:

-Upload a Seurat object file (.rds)
-Generate UMAP plot colored by metadata column
-Visualize expression of selected gene
-Download UMAP plot and gene expression plot
-Reset app state

## Installation
[App link](https://ndd0wk-sita0karan-patel.shinyapps.io/shiny2)

1. **Clone the repository:**
```
git clone https://github.com/your_username/your_app.git
cd your_app
```

2. **Run the Shiny app:**
```
```
library(shiny)
runApp("your_app_directory")
```
```

## How to Use This App

1. **Upload Seurat Object:**
   - Use the "Upload File" button to upload your Seurat object (in `.rds` format).
  
2. **Run Analysis:**
   - Press the "Run" button to initiate the analysis.
   - Explore UMAP plots and gene expression features in the generated tabs.
   - Press "Run" - this will generate the analysis
   - Select a metadata column to color the UMAP plot
   - Select a gene to visualize expression
   - Download plots using the download buttons

3. **Reset Inputs:**
   - To clear input files and remove all tabs, press the "Reset" button.

### Tab 1: Home Page
![Home Page](/images/home_page.png)

### Tab 2: scRNAseq Analyzer
![UMAP based on user selected metadata values](/images/scRNAseq_analyzer.png)
![Feauture plot of the user  seleted gene](/images/scRNAseq_analyzer.png)


## File Overview

The app consists of 3 main code files:

### app.R
- Defines the user interface (UI) and server logic.
- Implements the functionality for uploading files, running analysis, and creating tabs.
- Handles reactivity and renders outputs
- Buttons to run analysis, reset, and download plots

### global.R
- Loads required R packages
- Contains global functions and libraries.
- Defines functions like load_seurat_obj, create_metadata_UMAP, and create_feature_plot.
- Handles validation, UMAP plot creation, and feature plot generation.

### seurat.R
- Sample data - creates a Seurat object from raw Peripheral Blood Mononuclear Cells (PBMC) data from 10X Genomics
- Performs standard Seurat workflow (normalization, PCA, clustering, etc)
- Adds simulated metadata columns
- Saves Seurat object as a .rds file that can be used as input in the app.
Make sure to customize file paths, data sources, and additional processing steps according to your specific setup.

