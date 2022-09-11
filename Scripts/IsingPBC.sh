#!/bin/bash
#SBATCH --job-name=IsPBC2
#SBATCH --cpus-per-task=1      # ���-�� ���� �� ������ �������
#SBATCH --output=IsingPBC2.txt
#SBATCH --time=12:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

# ������������� ������ �����
srun /home/iipchelintsev/Ising/Ising.out 400 0.875 pbc &
srun /home/iipchelintsev/Ising/Ising.out 400 1.0 pbc &


# �������� ���������� ���� ���������
wait