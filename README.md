# ALPACA-paper

This repository contains instruction for reproducing figures from ALPACA manuscript


To test the ALPACA model, follow the instruction in https://github.com/McGranahanLab/ALPACA-model

To reproduce the figures follow these steps:
- [Using Conda](README.md#using-conda)
- [Using Docker image](README.md#using-docker-image)
- [Read this if you are using M1/M2 Apple Chip](README.md#apple-chips)

## Using Conda
1. Open the terminal and choose a directory, for example:
```
cd GitHub
```

2. If you don't have git installed, follow the instructions here: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

3. Clone the repository:
``` 
git clone https://github.com/McGranahanLab/ALPACA-paper.git
```

4. Navigate to the ALPACA-paper directory:
```
cd ALPACA-paper
```

5. Download the data package (~670 mb) from https://zenodo.org/records/15519765
```
curl -O https://zenodo.org/records/11060928/files/_assets.tar.gz
curl -O https://zenodo.org/records/11060928/files/output.tar.gz
```
or
```
wget https://zenodo.org/records/11060928/files/_assets.tar.gz
wget https://zenodo.org/records/11060928/files/output.tar.gz
```
You can also download the file via the internet browser and place the file in the ALPACA-paper directory

6. Extract the data to the ALPACA-paper directory
```
tar -xzvf _assets.tar.gz
tar -xzvf output.tar.gz
```
If you don't have `tar` command available, you can unpack the archive using 7-zip (https://www.7-zip.org)

7. Install conda if you don't already have it:

See instructions at: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

8. Create conda environment based on the supplied yaml file:

```
conda env create -f alpaca_figures.yml
```

9. If creating environment from the yml files does not work, create a new environment and add required packages and libraries OR see below for [Using Docker image](#using_docker_image)

```
conda create -n alpaca_figures python=3.8 r-essentials --channel conda-forge
conda activate alpaca_figures
# packages/libraries below are required to reproduce the figures:
pip install papermill
pip install jupyterlab
pip install seaborn
pip install plotly
pip install pandas
pip install kaleido
pip install networkx
conda install conda-forge::r-data.table
conda install conda-forge::r-dplyr
conda install conda-forge::r-ggpubr
conda install conda-forge::r-survminer
conda install conda-forge::r-survival
conda install conda-forge::r-optparse
conda install bioconda::bioconductor-genomicranges
conda install conda-forge::r-lmertest

# start R and install two remaining libraries not available via conda:
R # start R
install.packages("igraph")
install.packages("tidytext")
install.packages("ggridges")
q() # quit R


```
10. Activate the environment:
```
conda activate alpaca_figures
```

11. To reproduce the figures, navigate to project directory `cd ALPACA-paper` and execute `./run_all_figures.sh`. Make sure that `alpaca_figures` environment is active. Figures will be placed in `ALPACA-paper/figures` directory.

------
## Using Docker image:
- This procedure requires approximately 10GB of disk space to install containerised linux distribution
1. Install Docker and run Docker Desktop app: https://docs.docker.com/desktop/install/
2. Navigate to project directory `cd ALPACA-paper`
3. Build the image:


```
docker build -f Dockerfile -t alpaca_container .
```


4. Once the process is completed (~10 minutes) run the container:

```
docker run -it --name alpaca_test alpaca_container /bin/bash
```

5. You should now be within the containerised system and your terminal prompt should look similar to:

```
root@7e2a1eb2227a:/app#
```

6. To activate the conda virtual environment run the following commands:

```
conda init
source /root/.bashrc
conda activate alpaca
```

7. To recreate the figures run:

```
./run_all_figures.sh
```

8. To inspect the results of the model and figures, we need to export them from the container. First, exit the container by typing

```
exit
```

9. You should be back in the project root directory 'ALPACA-paper'. Copy the example output and figures with the following commands:

```
mkdir -p figures
docker cp alpaca_test:/app/figures .

mkdir -p output
docker cp alpaca_test:/app/output/example_cohort ./output
```

10. Stop and remove the image from your system with:
```
docker stop alpaca_test
docker rm alpaca_test
docker rmi alpaca_container
```

## Apple Chips
Docker images and conda environments might not work for Apple M1/M2 chips. In such case, install all the packages and libraries from the command line or R. Make sure to install at Python 3.9 and Gurobipy 11.
```
Python:
pip install pandas
pip install kneed
pip install gurobipy==11
# packages/libraries below are required to reproduce the figures:
pip install papermill
pip install jupyterlab
pip install seaborn
pip install plotly
pip install pandas
pip install kaleido
pip install networkx

R:
data.table
dplyr
ggpubr
survminer
survival
optparse
genomicranges
lmertest
igraph
tidytext
ggridges
```

------
This procedure has been tested in Linux (CentOS Linux 7, Linux 3.10.0-1160.62.1.el7.x86_64) and macOS (Sonoma 14.4.1) environments. For the full list of dependencies, please see the alpaca_figures.yml file. The test run, including downloads, environment creation and making the figures takes approximately 1-2 hrs on a standard laptop.
