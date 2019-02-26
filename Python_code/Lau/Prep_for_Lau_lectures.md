## Preparing for Lau lectures
You will need to use your Anaconda installation of Python and install a few other packages to run his code. If you have not already installed it, go here <https://www.anaconda.com/distribution/> and follow directions. If you want to learn more about Anaconda, check out this link <https://docs.anaconda.com/anaconda/user-guide/getting-started/>.

### Steps to generate a new conda environment and install required Python packages

1. Update the conda installation, create a new conda environment and activate it, and install other necessary packages. In shell (command line), type (or copy and paste) the commands below, one at a time. You will likely need to approve installation/updates along the way.

```
    conda update -n root conda
    conda create -n lau anaconda
    conda activate lau
    conda install -c anaconda jupyter
    pip install pydpc
```

Next, you will need to manually download some files to install one of the necessary packages (`scrublet`). This is from a Git repo on GitHub.

2. Download the scrublet git repo from GitHub <https://github.com/AllonKleinLab/scrublet> as a zip file (by default it usually gets saved into your Downloads folder) and uncompress (if not done automatically). You should see a new folder named `scrublet-master`.

3. Using the command line, change directories into the `scrublet-master` folder where you downloaded it onto your computer. Note, the `<>` symbols indicate you need to modify what is between them for the code to work on your computer. On my computer I would use `cd /Users/darren/Downloads/scrublet-master`

```{bash}
	cd <full-path-to>/scrublet-master
```

4. Once you are in the `scrublet-master` directory, type (or copy and paste) these commands:

```{bash}
	pip install -r requirements.txt
	pip install --upgrade .
```