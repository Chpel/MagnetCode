#!/bin/bash
#SBATCH --job-name=Dr_S_full
#SBATCH --output=Drunken_Sailor_colab_%a_k.txt
#SBATCH --time=2-00:00:00
#SBATCH --constraint="type_a"
#SBATCH --cpus-per-task=3
#SBATCH --array=1,2,5,10,50,100

module purge

# Загрузка необходимых переменных окружения
source ~/envs/colab/bin/activate

idx=$SLURM_ARRAY_TASK_ID
let "idx = idx * 100"
python3 -u /home/iipchelintsev/Magnet_Code/Random_Walk/Drunken_Sailor_v2.py $idx 10000 10000
