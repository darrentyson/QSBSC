# QSBSC
Quantitative Systems Biology of Single Cells class repository

We will add more here as the course proceeds.

## Installing R packages for Dr. Irish's lectures
If you have already cloned this repository onto your local computer, in R (or RStudio), you can type:
`source('<complete_path_to>/R_code/Irish_R_packages.r')` where `<complete_path_to>` is the file path to this repository on your own computer.

Alternatively, you can execute the code in R by accessing the file directly from GitHub using:
`source('https://github.com/darrentyson/QSBSS/raw/master/R_code/Irish_R_packages.r')`

## First homework assignment
Using your choice of how to interact with GitHub (Sourcetree, command line interface (CLI),
GitHub Desktop, etc.), you should clone this repository onto your local computer and perform 
the steps described in the commented lines at the bottom of 
`./R_code/Tyson/TestCodeForDataLoadingAndGit.r`

## Tip for keeping Git repo clutter-free
If you are using SourceTree you should go to Preferences and `Edit File` the Global Ignore list (a hidden file named `.gitignore_global` usually saved in your user directory that you could also edit manually.)

You should add things like this into the file:
```
*~
.DS_Store
.Rhistory
.Rapp.history
Iconr
*.pyc
```
So files matching these criteria do not get uploaded to GitHub when you commit (they will not be tracked).

## 2019-02-18: Added new flow cytometry data into `Data` directory
Please see the README file on GitHub https://github.com/darrentyson/QSBSC/tree/master/Data/CytobankData/experiment_46259_files
