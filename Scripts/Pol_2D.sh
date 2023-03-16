#!/bin/bash
#SBATCH --job-name=Pol_Tri_L
#SBATCH --output=Pol_Tri_Long_%A.txt
#SBATCH --time=7-00:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p forPolymer_Tri.txt)
echo $LINE
srun /home/iipchelintsev/Outs/$LINE
