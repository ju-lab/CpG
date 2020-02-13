#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/qsub_sdout/03_run_14.sh.sdout
cd /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision
python 03_reclassfication_add_individual_genotype.py 14