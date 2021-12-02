#!/bin/bash
#SBATCH --job-name=Is3D
#SBATCH --cpus-per-task=1      # ���-�� ���� �� ������ �������
#SBATCH --output=Is3D_Short.txt
#SBATCH --time=1-00:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3


for len in 100 300
do 
  srun /home/iipchelintsev/Ising3D.out $len 58 0 &
  srun /home/iipchelintsev/Ising3D.out $len 59 0 &
  for((i=61; i<69; i++))
  do
  srun /home/iipchelintsev/Ising3D.out $len $i 0 &
  done
done
wait