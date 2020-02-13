#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/qsub_sdout/01_annotation_long_run_by_chr_22.sh.sdout
cd /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision
python 01_annotation_long.py chr22