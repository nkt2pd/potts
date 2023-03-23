#!/bin/bash

#SBATCH --ntasks=1
#SBATCH -A cherngroup
#SBATCH --mem=32G
#SBATCH -p standard
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

./potts.exe
