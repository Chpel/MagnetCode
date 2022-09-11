#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # ���-�� ���� �� ������ �������
#SBATCH --output=Ising2D.txt
#SBATCH --time=1-00:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

# ������������� ������ �����

for len in 100 800 1600 3200
do 
  j=$(bc <<<"scale=3; 665 / 1000" )
  srun /home/iipchelintsev/TensorIn.out $len $j 1000000000 &
done


# �������� ���������� ���� ���������
wait