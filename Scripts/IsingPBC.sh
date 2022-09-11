#!/bin/bash
#SBATCH --job-name=IsPBC2
#SBATCH --cpus-per-task=1      # Кол-во ядер на каждый процесс
#SBATCH --output=IsingPBC2.txt
#SBATCH --time=12:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
module load gnu8 openmpi3

# Одновременный запуск задач
srun /home/iipchelintsev/Ising/Ising.out 400 0.875 pbc &
srun /home/iipchelintsev/Ising/Ising.out 400 1.0 pbc &


# Ожидание завершения всех процессов
wait