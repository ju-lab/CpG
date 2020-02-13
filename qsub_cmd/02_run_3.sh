#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/qsub_sdout/02_run_3.sh.sdout
cd /home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision
python 02_CGI_reclassificaiton_Anno.py 3