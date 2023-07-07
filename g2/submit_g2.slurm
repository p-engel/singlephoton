#!/bin/bash                                            
## Job Name                                            
#SBATCH --job-name=g2
## Allocation Definition                                      
#SBATCH --account=chem
#SBATCH --partition=chem
## Resources                                            
## Nodes                                              
#SBATCH --nodes=1                                         
## Tasks per node (28 is Slurm default) 
#SBATCH --ntasks-per-node=28 
## Walltime (days-HH:MM:SS)  
#SBATCH --time=3-07:00:00 
## Memory per node
#SBATCH --mem=240G                                        
##Output file                                           
#SBTACH --output output.out

module load anaconda3_5.3
python result.py

