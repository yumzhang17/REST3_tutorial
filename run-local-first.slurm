#!/bin/bash
#SBATCH --job-name="p53"        # job name
#SBATCH --nodes=2                # number of nodes
#SBATCH --gres=gpu:4             # number of gpus per node
#SBATCH -p fast                  # partition/queue
#SBATCH -n 64                   # number of cores
#SBATCH -o job%j.log   # without -e, combine stdout/err
#SBATCH -e job%j.err

now=$PWD
cd $now

set -x
date
at=`date +%s`

# gromacs
export PATH=/home/zgjia/Software/gromacs/gromacs2019_plumed252/bin:$PATH
export LD_LIBRARY_PATH=/home/zgjia/Software/gromacs/gromacs2019_plumed252/lib64:$LD_LIBRARY_PATH
# plumed
export PATH=/home/zgjia/Software/gromacs/plumed252/bin:$PATH
export LD_LIBRARY_PATH=/home/zgjia/Software/gromacs/plumed252/lib:$LD_LIBRARY_PATH
# mpi
module unload openmpi/gcc/64/1.10.1
export PATH=/home/zgjia/Software/openmpi/2.1.0/bin:$PATH
export LD_LIBRARY_PATH=/home/zgjia/Software/openmpi/2.1.0/lib:$LD_LIBRARY_PATH
# cuda
module load cuda/9.0.176

newdir=./

# make host file
hostfile=$now/job-host-$SLURM_JOB_ID
nodelist=`scontrol show hostname $SLURM_JOB_NODELIST`
echo "$nodelist" > $hostfile 
sed -i 's/$/ slots=32/g' $hostfile

# 8 replicas
# request 2 nodes => 4 replicas on each node
# for each replica, there is 1 mpi process, 8 threads used
# on one node, each replica use 1 of the 4 gpus
nsteps=500000000 # 1000 ns
mpirun -np 8 -npernode 4 -hostfile $hostfile gmx_mpi mdrun -plumed $now/plumed.dat -multidir $newdir/rep0 $newdir/rep1 $newdir/rep2 $newdir/rep3 $newdir/rep4 $newdir/rep5 $newdir/rep6 $newdir/rep7  -replex 5000 -nsteps $nsteps -hrex -s topol.tpr -deffnm prod

# write the time to a file
gmx_mpi check -f rep0/prod.cpt >& $now/check$$.log
frame1=`grep 'Last frame' $now/check$$.log | awk '{print $5}'`
rm -f $now/check$$.log
echo "$frame1" > $now/job-last-frame.dat

cd $now

date
bt=`date +%s`
duration1=`echo "($bt-$at)/60"   | bc -l`
duration2=`echo "($bt-$at)/3600" | bc -l`
echo "Run time is $duration1 minutes (i.e., $duration2 hours)"
perf=`echo "$nsteps*2/1000000/$duration2*24" | bc -l`
echo "Performance is $perf ns/day"

exit
