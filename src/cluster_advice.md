# Cluster advice



## Setting up a conda env:


module load miniconda3-4.5.4-gcc-5.4.0-hivczbz

#NOTE: to install a specific version of a package, do package_name=<version_number> e.g. iris=2.0.0
#create an environment thing:
conda create python=3.6.8 -p /home/lb584/software/conda_envs/iris_env --copy

#install extra modules inside the environemnt:
conda install -y -c conda-forge iris=2.4.0 -p /home/lb584/software/conda_envs/iris_env --copy