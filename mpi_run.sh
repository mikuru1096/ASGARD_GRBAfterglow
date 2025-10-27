#!/bin/bash
#SBATCH --account renjia_mcmc
#SBATCH -o run%j.out
#SBATCH -e error%j.out
#SBATCH -J MCMC_rj
#SBATCH -p normal
#SBATCH --exclusive
#SBATCH -t 72:00:00		#指定作业最大运行24小时
#SBATCH --mem=128G		#占用节点全部内存
#SBATCH --nodelist comput3
#SBATCH -N 1
#SBATCH -n 256

# running the command

mpirun -np 32 python pymultinest_demo.py
