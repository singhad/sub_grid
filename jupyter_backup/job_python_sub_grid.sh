#! /bin/bash
#SBATCH --job-name=analysis
#SBATCH --partition=teyssier
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 

source activate root

export XDG_RUNTIME_DIR=""
export PATH="/home/cluster/mkrets/anaconda2/bin:$PATH"
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$HOME/local/lib:$DYLD_LIBRARY_PATH
export PATH=$HOME/bin:$PATH

echo "python run"

srun python -uW ignore sub_grid.py

exit
