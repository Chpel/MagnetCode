#!/bin/bash
#SBATCH --job-name=Pol_Test
#SBATCH --output=Pol_Test_%A.txt
#SBATCH --time=01:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p forPolymer_Test.txt)
echo $LINE
srun /home/iipchelintsev/Magnet_Code/IsingUniversal_Shrinked/$LINE
