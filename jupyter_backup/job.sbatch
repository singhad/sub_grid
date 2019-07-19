#! /bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=ipython
#SBATCH --partition=teyssier
#SBATCH --output log_%J.txt
#SBATCH --time=23:45:00
#SBATCH --mem=12G
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1

cd;

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 

let ipnport=4169
echo ipnport=$ipnport
ipnip=$(hostname -i)
echo ipnip=$ipnip

echo "---------------------------------------------------------------------"
echo "For remote usage:"
echo "ssh -N -L $ipnport:$ipnip:$ipnport amasin@cluster.s3it.uzh.ch"
echo "---------------------------------------------------------------------"
echo ""

export XDG_RUNTIME_DIR=""

export PATH="/home/cluster/mkrets/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$HOME/local/lib:$DYLD_LIBRARY_PATH
export PATH=$HOME/bin:$PATH

echo "--- ipython run: ------------------------------------------------->>>"
jupyter notebook --no-browser --ip 0.0.0.0 --port=$ipnport
#jupyter lab --no-browser --ip 0.0.0.0 --port=$ipnport

