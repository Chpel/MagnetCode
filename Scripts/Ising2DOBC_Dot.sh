#!/bin/bash
#SBATCH --job-name=Ising2D
#SBATCH --cpus-per-task=1      # ���-�� ���� �� ������ �������
#SBATCH --output=Ising2D.txt
#SBATCH --time=1-00:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

# ������������� ������ �����

for len in 10 28 40 57
do 
  j=$(bc <<<"scale=4; 6673 / 10000" )
  r=$(bc <<<"scale=3; 473 / 1000" )
  srun /home/iipchelintsev/Ising2D.out $len $j obc $r &
done


# �������� ���������� ���� ���������
wait