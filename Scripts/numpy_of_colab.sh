#!/bin/bash
#SBATCH --job-name=n_o_c
#SBATCH --output=n_o_c.txt
#SBATCH --time=00:05:00

module purge

# Загрузка необходимых переменных окружения

source activate ~/envs/colab/bin/activate

python3 -u checker.py
