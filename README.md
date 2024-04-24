# ALPACA-paper

This repository contains: 
1) the ALPACA model, 
2) example script and data to test the model, 
3) outputs required to reproduce the figures.


To test the model and reproduce the figures follow these steps:

1. Open the terminal and choose a direcory, for example:
```
cd GitHub
```

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
wget https://zenodo.org/records/11060928/files/alpaca_data.tar.gz
```

6. Extract the data to the ALPACA-paper directory
```
tar -xzvf alpaca_data.tar.gz
```

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

9. Activate the environment:
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

9. Run the model for the example cohort:
```
python bin/MODEL/run_example.py
```

10. Expected output:
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

11. To reproduce the figures navigate to project directory `cd ALPACA-paper` and execute `./run_all_figures.sh`. Make sure that `alpaca_figures` or `alpaca` environemnts are active.
