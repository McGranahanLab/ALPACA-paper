# ALPACA-paper

This repository contains: 
1) the ALPACA model, 
2) example script and data to test the model, 
3) outputs required to reproduce the figures.


To test the model and reproduce the figures follow these steps:
- [Using Conda](#using_conda)
- [Using Docker image](#using_docker_image)

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

5. Download the data package (~670 mb) from https://zenodo.org/records/11060928/files/alpaca_data.tar.gz
```
curl -O https://zenodo.org/records/11060928/files/alpaca_data.tar.gz
```
or
```
wget https://zenodo.org/records/11060928/files/alpaca_data.tar.gz
```
You can also download the file via the internet browser and place the file in the ALPACA-paper directory

6. Extract the data to the ALPACA-paper directory
```
tar -xzvf alpaca_data.tar.gz
```
If you don't have `tar` command available, you can unpack the archive using 7-zip (https://www.7-zip.org)

7. Install conda if you don't already have it:

See instructions at: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

8. Create conda environment based on the supplied yaml files:

To just run the model:
```
conda env create -f alpaca_model.yml
```

To just recreate the figures:
```
conda env create -f alpaca_figures.yml
```

Or to do both:
```
conda env create -f alpaca.yml
```

9. If creating environment from the yml files does not work, create a new environment and add required packages and libraries OR see below for [Using Docker image](#using_docker_image)

```
conda create -n alpaca python=3.8 r-essentials
conda activate alpaca
pip install pandas #required for ALPACA model
pip install kneed #required for ALPACA model
pip install gurobipy==11 #required for ALPACA model
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
conda activate alpaca_model
```
or
```
conda activate alpaca_figures
```
or
```
conda activate alpaca
```

11. Run the model for the example cohort:
```
python bin/MODEL/run_example.py
```
The example run should take less than 5 minutes and once it is finished you should see:
```
Optimal solution found (tolerance 1.00e-04)
Best objective 0.000000000000e+00, best bound 0.000000000000e+00, gap 0.0000%

---------------------------------------------------------------------------
Multi-objectives: solved in 0.07 seconds (0.02 work units), solution count 9

Segment ALPACA_input_table_LTXSIM039_7_193245_2255974.csv done.
Creating combined output
Done

 _____ __    _____ _____ _____ _____
|  _  |  |  |  _  |  _  |     |  _  |
|     |  |__|   __|     |   --|     |
|__|__|_____|__|  |__|__|_____|__|__|
  /\⌒⌒⌒/\
  (⦿   ⦿)
  ( 'Y' )
   (   )
   (   )
   (   )
   (~ ~~~~~~~~~~)
   ( ~ ~~   ~~  )
   ( ~  ~ ~  ~  )
   (~  ~~~~~   ~)
    │ │     │ │
    │ │     │ │

```

12. Expected output:
Directory called 'output' will be created in the project directory and will contain the following files:
```
output
└── example_cohort
         ├── cohort_outputs
         │         └── combined.csv
         └── patient_outputs
                 ├── final_LTXSIM039.csv
                 └── final_LTXSIM001.csv
```
Each of these files contain clone specific copy number for each clone and each segment in the provided example cohort (inputs can be found in input directory).

13. To reproduce the figures, navigate to project directory `cd ALPACA-paper` and execute `./run_all_figures.sh`. Make sure that `alpaca_figures` or `alpaca` environments are active. Figures will be placed in `ALPACA-paper/figures` directory.

------
## Using Docker image:
- This procedure requires approximately 10GB of disk space to install containerised linux distribution
1. Install docker: https://docs.docker.com/desktop/install/
2. Navigate to project directory `cd ALPACA-paper`
3. Build the image:

```
docker build -t alpaca_container .
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

7. To run the model with the test cohort:

```
python bin/MODEL/run_example.py
```
See poin 11 above for the expexted output

8. To recreate the figures run:

```
./run_all_figures.sh
```

9. To inspect the results of the model and figures, we need to export them from the container. First, exit the container by typing

```
exit
```

10. You should be back in the project root directory 'ALPACA-paper'. Copy the example output and figures with the following commands:

```
mkdir -p figures
docker cp alpaca_test:/app/figures .

mkdir -p output
docker cp alpaca_test:/app/output/example_cohort ./output
```

11. Stop and remove the image from your system with:
```
docker stop alpaca_test
docker rm alpaca_test
docker rmi alpaca_container
```

------
This procedure has been tested in Linux (CentOS Linux 7, Linux 3.10.0-1160.62.1.el7.x86_64) and macOS (Sonoma 14.4.1) environments. For the full list of dependencies, please see the alpaca.yml file. The test run, including downloads, environment creation and making the figures takes approximately 1-2 hrs on a standard laptop.
