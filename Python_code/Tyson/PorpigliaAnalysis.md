#### Code to visualize Porpiglia et al. CyToF data (FR-FCM-ZY3K from FlowRepository.org)
Please use the `lau` conda environment by typing `conda activate lau` in the Anaconda Prompt or using the Anaconda Navigator app.

You will need to ensure that `cellxgene` is installed. You can check by typing `conda list` and searching for `cellxgene`.

If you need to install you can type `pip install cellxgene`. See https://pypi.org/project/cellxgene/ for more information.

The data has been processed first using the R code in `../../R_code/Porpiglia_paper/Porpiglia_analysis.r`

The data were then processed in `Python` using the code in `Porpiglia data cellxgene analysis.ipynb` and then using the following in the command line: `cellxgene prepare --layout umap ../../Data/Porpiglia/ZY3K.h5ad --output=../../Data/Porpiglia/ZY3K_umap.h5ad`. 

To visualize these data using `cellxgene` you should type the following in the command line (e.g. `Terminal`) 

```{/bin/bash}
cellxgene launch ../../Data/Porpiglia/ZY3K_umap.h5ad --open
```
You may need to change the path to the file based on whatever directory you are in at the command line. The path above is correct if you are in the same directory as this file.
