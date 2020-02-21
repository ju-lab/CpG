#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, os, gzip
chr_num = sys.argv[1]
#chr_num = '1'
output_file = file("signature_1000_chr%s.txt" % chr_num.replace("23","X"),"w")


if chr_num == '23' or chr_num == 'X':
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")
else:
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" %chr_num.replace("23","X"))   

info_line = info_file.readline().strip()
while info_line[0:1] == '#':
    info_line = info_file.readline().strip()

output_file.write('#CHROM\tPOS\tID\tREF\tALT\tAC\n')

while info_line:
    info_split = info_line.split('\t')
    info_chr = info_split[0]
    info_pos = int(info_split[1])
    info_ref = info_split[3]

    info_alt = info_split[4].split(',')
    for i in info_split[7].split(';'):
        if i.split('=')[0] == 'AC':
            info_ac = i.split('=')[1].split(',')      
    
    for j in info_alt:
        j_ac = int(info_ac[info_alt.index(j)])
        if  j_ac < 44:
            output_file.write("%s\t%s\t.\t%s\t%s\t%s\n" %(info_chr,info_pos,info_ref,j,j_ac))

    info_line = info_file.readline().strip()
output_file.close()


# In[ ]:




