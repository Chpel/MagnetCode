#!/bin/bash
for ((i=40 ; i <51 ; i+=1))
do
  for len in 100 300 400 500 1000
  do
  sbatch --time=5-0:0 --wrap="srun ./Ising3D.out $len $i 0"
  done
done