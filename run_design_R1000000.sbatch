#!/bin/bash
#SBATCH -N 1 -n 100 -t 24:00:00 -C rome -p ccb --export=ALL
#SBATCH --job-name=Design_R1000000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qzhu@flatironinstitute.org

source /mnt/home/vmulligan/load_my_modules.sh
source /mnt/home/vmulligan/masala_workingcopy/set_up_masala.sh

mpirun -np 25 /mnt/home/vmulligan/rosetta_git_workingcopy/Rosetta/main/source/bin/rosetta_scripts.masalampiserialization.linuxgccrelease -masala_plugins ~vmulligan/masala_workingcopy/standard_masala_plugins/ -total_threads 4 -masala_total_threads 4 @design_R1000000.flags
