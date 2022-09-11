#!/bin/bash
#SBATCH --job-name=Sp_simple
#SBATCH --output=Sp_simple.txt
#SBATCH --time=1-00:00:00
#SBATCH --constraint="type_a"

module purge

# Загрузка необходимых переменных окружения
source ~/envs/colab/bin/activate

python3 -u /home/iipchelintsev/Magnet_Code/Random_Walk/Spitser_Simpler.py 10000 1000000000
