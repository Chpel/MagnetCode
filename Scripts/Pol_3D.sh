#!/bin/bash
#SBATCH --job-name=Pol_3D
#SBATCH --output=Pol_3D_%A.txt
#SBATCH --time=7-00:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
module load gnu8 openmpi3

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p forPolymer_3D.txt)
echo $LINE
srun /home/iipchelintsev/Outs/$LINE

