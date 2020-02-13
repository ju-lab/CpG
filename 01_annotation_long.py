#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#add thousand genome and PCAWG data including individual genotypes

import sys, os, gzip
target_chr = sys.argv[1] #ex> chr1, chrX, chrY

input_file = file("revision_CpG_CtoT.txt")
output_file = file("revision_CpG_CtoT_PCAWG_%s_individual.txt" % target_chr,"w")
output_file1 = file("revision_CpG_CtoT_PCAWG_%s.txt" % target_chr,"w")
input_line = input_file.readline().strip()

if target_chr == 'chrX':
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")
elif target_chr == 'chrY':
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz")     
else:
    info_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/1000Genomes/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % target_chr)   

pcawg_file = gzip.open("/home/users/jhyouk/89_backup_Workstation/jhyouk/CpG_revision/PCAWG-8/pcawg8.snps.indels.svs.phased.icgc.v2.vcf.gz")

info_line = info_file.readline().strip()
pcawg_line = pcawg_file.readline().strip()
while pcawg_line[0:2] == '##':
    pcawg_line = pcawg_file.readline().strip()
while info_line[0:2] == '##':
    info_line = info_file.readline().strip()
    
info_chr = 1

pcawg_split = pcawg_line.split('\t')
info_split = info_line.split('\t')

output_file.write(input_line + '\tPCAWG\tThousand\tTotal_Af;AC;AN\tRace_thousand_EAS;AMR;AFR;EUR;SAS\tGT_PCAWG;Thousand\t' + '\t'.join(pcawg_split[9:]) + '\t'+  '\t'.join(info_split[9:]) +'\n')
output_file1.write(input_line + '\tPCAWG\tThousand\tTotal_Af;AC;AN\tRace_thousand_EAS;AMR;AFR;EUR;SAS\tGT_PCAWG;Thousand' +'\n')

input_line = input_file.readline().strip()
info_line = info_file.readline().strip()
pcawg_line = pcawg_file.readline().strip()

p_done = 'no'


while input_line:
    input_split = input_line.split('\t')
    if input_split[0] == target_chr:
        break
    else:
        input_line = input_file.readline().strip()
    
    
while pcawg_line:
    p_split = pcawg_line.split('\t')
    if 'chr' + p_split[0] == target_chr:
        break
    else:
        pcawg_line = pcawg_file.readline().strip()
        
while info_line:
    info_split = info_line.split('\t')
    if 'chr' + info_split[0] == target_chr:
        break
    else:
        info_line = info_file.readline().strip()



while input_line:
    input_split = input_line.split('\t')
    if input_split[0] != target_chr:
        break
        
    input_chr = int(input_split[0].replace("chr","").replace("X","23").replace("Y","24"))
    input_pos = int(input_split[1])
    input_ref = input_split[3]
    input_alt = input_split[4]
    
    p_AC = 0
    p_AN = 0
    
    t_info = '.'
    race_info = '.'
    
    gt_pcawg = 0
    gt_thousand = 0   
    
    
    p_indi = []
    for i in range(0,1823):
        p_indi.append('.')
    t_indi = []
    for i in range(0,2504):
        t_indi.append('.')    
        
    #PCAWG annotation
    if p_done == 'no':
        p_info = '.'
        if pcawg_line == '':
            p_done = 'yes'
            p_info = '.;.;.'              
            
        else:

            p_split = pcawg_line.split('\t')

            if 'chr' + p_split[0] != target_chr:
                p_done = 'yes'
                p_info = '.;.;.'            

            p_chr = int(p_split[0].replace("chr","").replace("X","23").replace("Y","24"))
            p_pos = int(p_split[1])
            #dec_target = 'none'
            if p_chr < input_chr:
                pcawg_line = pcawg_file.readline().strip()
                continue
            elif p_chr == input_chr:
                if p_pos < input_pos:
                    pcawg_line = pcawg_file.readline().strip()
                    continue
                elif p_pos > input_pos:
                    p_done = 'yes'
                    p_info = '.;.;.'
                    #gt_pcawg = 0
                else:
                    p_ref = p_split[3]
                    p_alt = p_split[4].split(',')

                    p_target = 0
                    dec_target = '.'
                    if p_ref == 'C':
                        for p in p_alt:
                            if p == 'T':
                                dec_target = 'bingo'
                                break
                            else:
                                p_target+=1
                    elif p_ref == 'G':
                        for p in p_alt:
                            if p == 'A':
                                dec_target = 'bingo'
                                break
                            else:
                                p_target+=1
                    else:
                        'blank'

                    if dec_target == 'bingo':
                        gt_pcawg = p_target + 1
                        p_num = p_split[7].split(';')
                        for i in p_num:
                            if i[0:3] == 'AN=':
                                p_AN = float(i.split('=')[1])
                            elif i[0:3] == 'AC=':
                                temp_AC = i.split('=')[1].split(',')
                                p_AC = float(temp_AC[p_target])
                        p_info = '%s;%s;%s' %(round(p_AC/p_AN,5),p_AC,p_AN)
                        p_indi = p_split[9:]
                        #nn+=1
                    else:
                        p_info = '.;.;.'


                    p_done = 'yes'
                    pcawg_line = pcawg_file.readline().strip()

            else:
                p_done = 'yes'
                p_info = '.;.;.'
            
    else: #p_done == 'yes
        'blank'
    
    #thousand genome annotation
    info_AC = 0
    info_AN = 0
    info_race = '.'
    af_race = ['0','0','0','0','0']
    info_target = 'none'
    
    if info_line == '':
        t_info = '.;.;.'
        info_race = '.;.;.;.;.'        
    else:
        info_split = info_line.split('\t')        
        info_chr = int(info_split[0].replace("chr","").replace("X","23").replace("Y","24"))
        info_pos = int(info_split[1])

        if info_chr < input_chr:
            info_line = info_file.readline().strip()
            continue
        elif info_chr == input_chr:
            if info_pos < input_pos:
                info_line = info_file.readline().strip()
                continue            
            elif info_pos > input_pos:
                t_info = '.;.;.'
                info_race = '.;.;.;.;.'
            else:
                info_ref = info_split[3]
                info_alt = info_split[4].split(',')

                t_target =0
                t_decision = '.'
                if info_ref == 'C':
                    for t in info_alt:
                        if t == 'T':
                            t_decision = 'bingo'
                            break
                        else:
                            t_target+=1
                elif info_ref == 'G':
                    for t in info_alt:
                        if t == 'A':
                            t_decision = 'bingo'
                            break
                        else:
                            t_target +=1
                else:
                    'blank'

                if t_decision == 'bingo':
                    gt_thousand = t_target + 1
                    t_num = info_split[7].split(';')
                    for i in t_num:
                        #print input_line
                        #print i
                        if i[0:3] == 'AN=':
                            info_AN = float(i.split('=')[1])
                        elif i[0:3] == 'AC=':
                            temp_AC = i.split('=')[1].split(',')
                            info_AC = float(temp_AC[t_target])
                        elif i[0:3] == 'EAS':
                            temp= i.split('=')[1].split(',')
                            af_race[0] = str(temp[t_target])                        
                        elif i[0:3] == 'AMR':
                            temp= i.split('=')[1].split(',')
                            af_race[1] = str(temp[t_target])                          
                        elif i[0:3] == 'AFR':
                            temp= i.split('=')[1].split(',')
                            af_race[2] = str(temp[t_target])                          
                        elif i[0:3] == 'EUR':
                            temp= i.split('=')[1].split(',')
                            af_race[3] = str(temp[t_target])  
                        elif i[0:3] == 'SAS':
                            temp= i.split('=')[1].split(',')
                            af_race[4] = str(temp[t_target])  

                    t_info = '%s;%s;%s' %(round(info_AC/info_AN,5),info_AC,info_AN)   
                    info_race = ';'.join(af_race)
                    #print '\t'.join(info_split[9:15])
                    t_indi = info_split[9:]

                else:
                    t_info = '.;.;.'
                    info_race = '.;.;.;.;.'

        else:
            t_info = '.;.;.'
            info_race = '.;.;.;.;.'
    
    total_info = '.'
    total_AC = float(p_AC + info_AC)
    total_AN = float(p_AN + info_AN)
    total_af = 0
    if total_AN == 0:
        total_info = '.;.;.'
    else:
        total_af = round(float(total_AC/total_AN),5)
        total_info = str(total_af) + ';' + str(total_AC) + ';' + str(total_AN)
    
    
    output_file.write(input_line + '\t' + p_info +'\t' + t_info + '\t' + total_info + '\t' + info_race + '\t' + str(gt_pcawg) + ';' + str(gt_thousand) + '\t' + '\t'.join(p_indi) + '\t' + '\t'.join(t_indi) + '\n')
    output_file1.write(input_line + '\t' + p_info +'\t' + t_info + '\t' + total_info + '\t' + info_race + '\t' + str(gt_pcawg) + ';' + str(gt_thousand) + '\n')
    p_done = 'no'
    input_line = input_file.readline().strip()

    #nn+=1
    #if nn>1:
        #break

output_file.close()
output_file1.close()
info_file.close()
pcawg_file.close()
input_file.close()

print 'THE END'


