# Cluster advice for python users:

## Setting up a conda env:

Load anaconda
```bash
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
```

Create an anaconda virtual environment:

```bash
conda create python=3.9.0 -p /home/nt409/software/conda_envs/hrhr_env --copy
```
NB: to install a specific version of a package, do package_name=<version_number> e.g. iris=2.0.0

Install any python modules required inside the environment:

```bash
conda install -y -c conda-forge iris=2.4.0 -p /home/nt409/software/conda_envs/hrhr_env --copy
```

## Using a conda env:

To activate:
```bash
source activate /home/nt409/software/conda_envs/hrhr_env
```
To deactivate:
- open a new login, or

```bash
conda deactivate
```



## Script for job:

Set up with boilerplate from wiki.

Removed:
```bash
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module load use.own                        # This line loads the own module list
```

Added:
```bash
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate /home/nt409/software/conda_envs/hrhr_env
```


## To test a job:

```bash
bash foldername/myjobname
```

or just run individual python file locally:
```bash
python -m foldername.myscriptname PARAM1
```

## To submit a job:

```bash
sbatch foldername/myjobname
```

e.g.
```bash
sbatch param_scan/scan.submit
```


To check progress:

```bash
squeue -u nt409
```

Est start time:

```bash
squeue -u nt409 --start
```

