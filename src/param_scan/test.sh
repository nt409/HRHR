module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
module load use.own                        # This line loads the own module list
module load /rds/project/cag1/rds-cag1-general/epidem-modules/epidem.modules         # Loads Epidemiology group module list
module load miniconda3/4.9.2

# Conda set up
# >>> conda initialize >>>
# Contents within this block are managed by 'conda init' !!
__conda_setup="$('/rds/project/cag1/rds-cag1-general/epidem-programs/miniconda3/4.9.2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/rds/project/cag1/rds-cag1-general/epidem-programs/miniconda3/4.9.2/etc/profile.d/conda.sh" ]; then
        . "/rds/project/cag1/rds-cag1-general/epidem-programs/miniconda3/4.9.2/etc/profile.d/conda.sh"
    else
        export PATH="/rds/project/cag1/rds-cag1-general/epidem-programs/miniconda3/4.9.2/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate /home/nt409/software/conda_envs/hrhr_env

python -m param_scan.run.r4_re_run_failures_cluster 0