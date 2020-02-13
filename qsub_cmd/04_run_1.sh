#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/qsub_sdout/04_run_1.sh.sdout
cd /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision
python 04_union_indivial.py &> 04_run.out