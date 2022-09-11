#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # Кол-во ядер на каждый процесс
#SBATCH --output=Ising2D.txt
#SBATCH --time=1-00:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
module load gnu8 openmpi3

# Одновременный запуск задач

for len in 100 800 1600 3200
do 
  j=$(bc <<<"scale=3; 665 / 1000" )
  srun /home/iipchelintsev/TensorIn.out $len $j 1000000000 &
done


# Ожидание завершения всех процессов
wait