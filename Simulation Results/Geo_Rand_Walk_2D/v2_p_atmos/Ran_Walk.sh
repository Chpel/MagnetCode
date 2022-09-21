#!/bin/bash
#SBATCH --job-name=Dr_S_full
#SBATCH --output=Drunken_Sailor_%A.txt
#SBATCH --time=2-00:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
#SBATCH --cpus-per-task=5
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
source ~/envs/colab/bin/activate

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p forRanWalk.txt)
echo $LINE
python3 -u /home/iipchelintsev/Magnet_Code/Random_Walk/$LINE