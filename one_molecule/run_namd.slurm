#!/bin/bash
#SBATCH --job-name="n_m"  #  Please attempt to give the job a descriptive name
#SBATCH --partition=gpu,pi.cisborn,cenvalarc.gpu
#SBATCH --nodes=1                    # must be 1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1                    # max 4.5 times the number of tasks per node
#SBATCH --mem=16G
#SBATCH --time=02:15:00
#SBATCH --export=ALL



slurm_info_out(){
echo "=================================== SLURM JOB ==================================="
date
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "=================================== SLURM JOB ==================================="
echo "--- SLURM job-script output ---"
}

slurm_startjob () {
slurm_info_out


# Load Terachem Burrata
#export PATH="/data/cmyers7/terachem_2023:$PATH"

# Python and CUDA stuff
#module purge
#module load anaconda3
#CUDA=$(module -t spider cuda 2>/dev/null | grep -Eo 'cuda/[0-9.]+' | sort -V | tail -1 || true)
#if [[ -n $CUDA ]]; then
#  echo "Loading $CUDA"
#  module load "$CUDA"
#else
#  echo "No standalone CUDA module; relying on Anaconda’s cudatoolkit."
#fi


# Load Terachem module load cuda/11.4.0
module load anaconda3 >> env.out
source ~/.bashrc >> env.out 
conda activate pysces  >> env.out
#module load cuda/12.3.0

export TeraChem=/home/cmyers7/code/terachem/build
export PATH=$TeraChem/bin:$PATH
export PATH=/home/cmyers7/code/terachem/build/bin/:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/cmyers7/code/terachem/build/lib/
export TeraChem=/home/cmyers7/code/terachem/build
ldd /home/cmyers7/code/terachem/build/bin/terachem >> env.out
which terachem >>  env.out
#echo $(pwd) >> path.out

sleep 2

echo "starting jobs" >> status.out
for i in {1..40}
do    
   bash setup.sh   
   cd traj_1
   cd step
   echo $i >> ../../status.out
   pwd >> ../../status.out 
   cp step.xyz coordinate.xyz
   terachem input.in > output.out
#   sleep 200
   cd ..
   cd ..
#   cp step.xyz coordinate.xyz 
#   terachem input.in > output.out 
#   sleep 25
#   cd ..
    pwd >> status.out 
    bash extract_traj_data.sh >> a.out
    bash time_data.sh >> a.out
#    python propagate_traj.py  >> status.out 
    echo "before propagate" >> status.out
    python propagate_traj.py >> status.out 2>&1
    echo "after propagate, exit code = $?" >> status.out
done 

echo "Done" >> status.out


echo Job Done
}

slurm_startjob

