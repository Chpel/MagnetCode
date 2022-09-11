#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # ���-�� ���� �� ������ �������
#SBATCH --output=Ising2D.txt
#SBATCH --time=12:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

# ������������� ������ �����

i=$(bc <<<"scale=3; 463 / 1000" )
srun /home/iipchelintsev/Ising2D.out 1000 $i &

# �������� ���������� ���� ���������
wait