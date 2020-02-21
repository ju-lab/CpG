#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys, os, gzip
chr_num = sys.argv[1]
#chr_num = '1'
output_file = file("Cosmic_signature_extraction_chr%s.txt" % chr_num.replace("23","X"),"w")


if chr_num == '23' or chr_num == 'X':
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")
else:
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" %chr_num.replace("23","X"))   

pcawg_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/PCAWG-8/pcawg8.snps.indels.svs.phased.icgc.v2.vcf.gz")
info_line = info_file.readline().strip()
pcawg_line = pcawg_file.readline().strip()
while pcawg_line[0:1] == '#':
    pcawg_line = pcawg_file.readline().strip()
while info_line[0:1] == '#':
    info_line = info_file.readline().strip()

print 'position0'

while pcawg_line:
    p_split = pcawg_line.split('\t')
    if p_split[0] == chr_num.replace("23","X"):
        break
        
    pcawg_line = pcawg_file.readline().strip()    

print 'position1'
output_file.write('#CHROM\tPOS\tID\tREF\tALT\tAC\tAF\n')

info_dic = {}
m=0
while info_line:
    info_split = info_line.split('\t')
    info_chr = info_split[0]
    info_pos = int(info_split[1])
    info_ref = info_split[3]

    info_alt = info_split[4].split(',')
    for i in info_split[7].split(';'):
        if i.split('=')[0] == 'AC':
            info_ac = i.split('=')[1].split(',')      
    
    info_dic[info_pos] = [info_ref,info_alt,info_ac]
    info_line = info_file.readline().strip()
    m+=1

print 'position2'
m=0
p_dic = {}
while pcawg_line:
    p_split = pcawg_line.split('\t')
    p_chr = p_split[0]
    if p_chr != chr_num.replace("23","X"):
        break
    p_pos = int(p_split[1])
    p_ref = p_split[3]

    p_alt = p_split[4].split(',')
    for j in p_split[7].split(';'):
        if j.split('=')[0] == 'AC':
            p_ac = j.split('=')[1].split(',')     
        
    p_dic[p_pos] = [p_ref,p_alt,p_ac]
    pcawg_line = pcawg_file.readline().strip() 
    m+=1

        
print len(info_dic)
print len(p_dic)
input_min = min(min(info_dic.keys()),min(p_dic.keys()))
input_max = max(max(info_dic.keys()),max(p_dic.keys()))

ref_list = ['A','C','T','G']

for input_pos in range(input_min,input_max+1):
    info_ref = '.'
    info_alt = '.'
    info_ac = 0
    p_ref = '.'
    p_alt = '.'
    p_ac = 0
    
    if input_pos in info_dic.keys():
        if info_dic[input_pos][0] in ref_list:
            info_ref = info_dic[input_pos][0]
            info_alt = info_dic[input_pos][1]
            info_ac = info_dic[input_pos][2]
    
    if input_pos in p_dic.keys():
        if p_dic[input_pos][0] in ref_list:
            p_ref = p_dic[input_pos][0]
            p_alt = p_dic[input_pos][1]
            p_ac = p_dic[input_pos][2]
    true_ref = '.'
    if info_ac == 0 and p_ac == 0:
        continue
    elif info_ac > 0:
        true_ref = info_ref
    elif p_ac > 0:
        true_ref = p_ref
        
    temp_list = [0,0,0,0]
        
    for ii in info_alt:
        if ii in ref_list:
            temp_list[ref_list.index(ii)] += int(info_ac[info_alt.index(ii)])
    for jj in p_alt:
        if jj in ref_list:
            try:
                temp_list[ref_list.index(jj)] += int(p_ac[p_alt.index(jj)])    
            except:
                print input_pos
                print p_alt
                print p_ac
                sys.exit(1)
    
    for k in range(0,4):
        if temp_list[k] > 0:
            if float(temp_list[k])/4327.0 <=0.01:
                output_file.write("%s\t%s\t.\t%s\t%s\t%s\t%s\n" % (chr_num.replace("23","X"),input_pos,true_ref,ref_list[k],temp_list[k],float(temp_list[k])/4327.0,))


        
output_file.close()             
print 'THE END'     

