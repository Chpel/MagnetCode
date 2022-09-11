#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # Кол-во ядер на каждый процесс
#SBATCH --output=Ising2D.txt
#SBATCH --time=12:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
module load gnu8 openmpi3

# Одновременный запуск задач

i=$(bc <<<"scale=3; 463 / 1000" )
srun /home/iipchelintsev/Ising2D.out 1000 $i &

# Ожидание завершения всех процессов
wait