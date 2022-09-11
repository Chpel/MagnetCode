#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # Кол-во ядер на каждый процесс
#SBATCH --output=Ising2D.txt
#SBATCH --time=1-00:00:00                  # Ограничение времени выполнения задачи (дни-часы:мин:сек)
module purge                             # Очистка переменных окружения

# Загрузка необходимых переменных окружения
module load gnu8 openmpi3

# Одновременный запуск задач

for len in 10 28 40 57
do 
  j=$(bc <<<"scale=4; 6673 / 10000" )
  r=$(bc <<<"scale=3; 473 / 1000" )
  srun /home/iipchelintsev/Ising2D.out $len $j obc $r &
done


# Ожидание завершения всех процессов
wait