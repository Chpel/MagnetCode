#!/bin/bash
#SBATCH --job-name=Pol_4D
#SBATCH --output=Pol_4D_%A.txt
#SBATCH --time=3-00:00:00                  # ����������� ������� ���������� ������ (���-����:���:���)
module purge                             # ������� ���������� ���������

# �������� ����������� ���������� ���������
module load gnu8 openmpi3

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p forPolymer_4D.txt)
echo $LINE
srun /home/iipchelintsev/Outs/$LINE
